#!/bin/bash

# Default query (if user does not provide one)
SEARCH_QUERY="Human Microbiome Project AND amplicon[Strategy] AND Illumina[Platform]"

# Function to show usage
usage() {
    echo "Usage: $0 [-q \"<NCBI Search Query>\"]"
    echo "Example: $0 -q \"Human Microbiome Project AND amplicon[Strategy] AND Illumina[Platform]\""
    exit 1
}

# Parse command-line options
while getopts "q:" opt; do
    case $opt in
        q) SEARCH_QUERY="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure R is installed via Conda, install if missing
if ! command -v R &> /dev/null; then
    echo "R is not installed. Installing R using Conda..."
    conda install -y r-base
fi

ESCAPED_QUERY=$(echo "$SEARCH_QUERY" | sed "s/\"/\\\\\"/g")

# Create the R script dynamically with the user-provided query
cat <<EOF > hmp_16s_query.R
# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load required libraries, install if missing
required_packages <- c("rentrez", "dplyr", "stringr", "purrr", "xml2", "tidyr", "readr", "foreach", "doParallel", "future", "furrr", "parallelly")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Detect available CPU cores, respecting system limits
available_cores <- parallelly::availableCores()
message("Available CPU cores: ", available_cores)

# If only 1 core is available, run sequentially
if (available_cores <= 1) {
  message("Only 1 core available. Running sequentially...")
  plan(sequential)  # Use sequential processing
} else {
  # Use fewer cores than available to avoid overloading the system
  num_cores <- min(available_cores - 1, 4)  # Use at most 4 cores
  message("Using ", num_cores, " cores for parallel processing...")
  plan(multisession, workers = num_cores)  # Use future's multisession for parallel processing
}

batch_size <- 10000  # Adjustable
delay_time <- 3  # Increased delay to avoid API limits
max_retries <- 5  # Maximum number of retries for failed batches
search_query <- "$ESCAPED_QUERY"

search_results <- entrez_search(
  db = "sra",
  term = search_query,
  retmax = 0,
  use_history = TRUE
)

total_records <- search_results\$count
message("Total Records Found: ", total_records)

# Function to fetch and parse a batch of records with retry logic
fetch_batch <- function(start, batch_size, search_results, db_name) {
  retries <- 0
  while (retries < max_retries) {
    tryCatch({
      metadata_xml <- entrez_fetch(db = db_name, web_history = search_results\$web_history, rettype = "xml", retstart = start, retmax = batch_size)
      parsed_metadata <- read_xml(metadata_xml)
      
      # Extract metadata fields dynamically
      all_nodes <- xml_find_all(parsed_metadata, "//*")  # Get all XML nodes
      field_names <- unique(xml_name(all_nodes))  # Extract unique field names
      
      # Store extracted data
      batch_data <- list()
      for (field in field_names) {
        batch_data[[field]] <- xml_find_all(parsed_metadata, paste0("//", field)) %>% xml_text()
      }
      
      # Standardize all fields to the same length
      max_length <- max(sapply(batch_data, length))
      for (field in names(batch_data)) {
        batch_data[[field]] <- c(batch_data[[field]], rep(NA, max_length - length(batch_data[[field]])))
      }
      
      # Convert batch data to dataframe
      batch_metadata <- as_tibble(batch_data)
      return(batch_metadata)
    }, error = function(e) {
      retries <- retries + 1
      cat("Retry", retries, "for batch starting at", start, "due to error:", e\$message, "\\n")
      Sys.sleep(delay_time * retries)  # Increase delay with each retry
    })
  }
  warning("Failed to fetch batch starting at", start, "after", max_retries, "retries.")
  return(NULL)
}

# Batch Processing (parallel or sequential depending on available cores)
batch_indices <- seq(0, total_records, by = batch_size)

if (available_cores <= 1) {
  # Sequential processing
  all_metadata <- map_dfr(batch_indices, function(start) {
    cat("Fetching records", start, "to", min(start + batch_size, total_records), "\\n")
    result <- fetch_batch(start, batch_size, search_results, "sra")
    if (is.null(result)) {
      cat("Failed to fetch batch starting at", start, "after", max_retries, "retries. Skipping...\\n")
      return(NULL)
    }
    return(result)
  })
} else {
  # Parallel processing with future
  all_metadata <- future_map_dfr(batch_indices, function(start) {
    cat("Fetching records", start, "to", min(start + batch_size, total_records), "\\n")
    result <- fetch_batch(start, batch_size, search_results, "sra")
    if (is.null(result)) {
      cat("Failed to fetch batch starting at", start, "after", max_retries, "retries. Skipping...\\n")
      return(NULL)
    }
    return(result)
  }, .options = furrr_options(seed = TRUE))
}

# Validate total rows
if (nrow(all_metadata) != total_records) {
  warning("Mismatch in total records. Expected:", total_records, "Got:", nrow(all_metadata))
} else {
  cat("All records fetched successfully. Total rows:", nrow(all_metadata), "\\n")
}

# Save raw data
write_csv(all_metadata, "raw_sra_metadata.csv")

cat("Raw metadata saved successfully: raw_sra_metadata.csv\\n")
EOF

# Run the R script
Rscript hmp_16s_query.R

echo "Script execution completed."
