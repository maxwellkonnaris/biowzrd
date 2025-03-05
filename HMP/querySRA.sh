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

message("Starting R script execution...")
flush.console()

# Define required packages
required_packages <- c("rentrez", "dplyr", "stringr", "purrr", "xml2", 
                       "tidyr", "readr", "foreach", "doParallel", 
                       "future", "furrr", "parallelly")

# Identify missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

# Install missing packages if any
if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    flush.console()
    install.packages(missing_packages, dependencies = TRUE)
}

# Load all packages and check if any failed
loaded_packages <- sapply(required_packages, function(pkg) require(pkg, character.only = TRUE))

# Ensure furrr is properly loaded
if (!loaded_packages["furrr"]) {
    message("furrr package not found. Installing now...")
    flush.console()
    install.packages("furrr", dependencies = TRUE)
    library(furrr)
}

# Warn the user if any package failed to load
if (any(!loaded_packages)) {
    warning("The following packages failed to load: ", 
            paste(required_packages[!loaded_packages], collapse = ", "))
    flush.console()
}

# Detect available CPU cores, respecting system limits
available_cores <- parallelly::availableCores()
message("Available CPU cores: ", available_cores)
flush.console()

# If only 1 core is available, run sequentially
if (available_cores <= 1) {
  message("Only 1 core available. Running sequentially...")
  flush.console()
  plan(sequential)  # Use sequential processing
  batch_size <- 10000  # Larger batch size for sequential processing
} else {
  # Use fewer cores than available to avoid overloading the system
  num_cores <- min(available_cores - 1, 4)  # Use at most 4 cores
  message("Using ", num_cores, " cores for parallel processing...")
  flush.console()
  plan(multisession, workers = num_cores)  # Use future's multisession for parallel processing
  batch_size <- 1000  # Smaller batch size for parallel processing
}

delay_time <- 3  # Increased delay to avoid API limits
max_retries <- 5  # Maximum number of retries for failed batches
search_query <- "$ESCAPED_QUERY"

message("Executing search query: ", search_query)
flush.console()

search_results <- entrez_search(
  db = "sra",
  term = search_query,
  retmax = 0,
  use_history = TRUE
)

total_records <- search_results\$count
message("Total Records Found: ", total_records)
flush.console()

# Function to fetch and parse a batch of records with retry logic
fetch_batch <- function(start, batch_size, search_results, db_name) {
  message("Attempting to fetch batch: ", start, " - ", min(start + batch_size, total_records))
  flush.console()
  retries <- 0
  while (retries < max_retries) {
    tryCatch({
      metadata_xml <- entrez_fetch(db = db_name, web_history = search_results$web_history, rettype = "xml", retstart = start, retmax = batch_size)
      parsed_metadata <- read_xml(metadata_xml)
      
      if (is.null(parsed_metadata)) {
        message("Parsed metadata is NULL for batch: ", start)
        flush.console()
        return(NULL)
      }
      
      message("Successfully fetched batch: ", start)
      flush.console()

      # Extract metadata fields dynamically
      all_nodes <- xml_find_all(parsed_metadata, "//*")  # Get all XML nodes
      field_names <- unique(xml_name(all_nodes))  # Extract unique field names

      batch_data <- list()
      for (field in field_names) {
        batch_data[[field]] <- xml_find_all(parsed_metadata, paste0("//", field)) %>% xml_text()
      }

      # Standardize all fields to the same length
      max_length <- max(sapply(batch_data, length))
      for (field in names(batch_data)) {
        batch_data[[field]] <- c(batch_data[[field]], rep(NA, max_length - length(batch_data[[field]])))
      }

      batch_metadata <- as_tibble(batch_data)
      return(batch_metadata)
    }, error = function(e) {
      retries <- retries + 1
      message("Retry ", retries, " for batch at ", start, " due to error: ", e$message)
      flush.console()
      Sys.sleep(delay_time * retries)  # Increase delay with each retry
    })
  }
  warning("Failed to fetch batch starting at ", start, " after ", max_retries, " retries.")
  flush.console()
  return(NULL)
}

batch_indices <- seq(0, total_records - 1, by = batch_size)

if (available_cores <= 1) {
  all_metadata <- map_dfr(batch_indices, function(start) {
    message("Fetching records ", start, " to ", min(start + batch_size, total_records))
    flush.console()
    result <- fetch_batch(start, batch_size, search_results, "sra")
    return(result)
  })
} else {
  all_metadata <- future_map_dfr(batch_indices, function(start) {
    message("Fetching records ", start, " to ", min(start + batch_size, total_records))
    flush.console()
    result <- fetch_batch(start, batch_size, search_results, "sra")
    return(result)
  }, .options = furrr_options(seed = TRUE))
}

if (nrow(all_metadata) != total_records) {
  warning("Mismatch in total records. Expected:", total_records, "Got:", nrow(all_metadata))
  flush.console()
} else {
  message("All records fetched successfully. Total rows: ", nrow(all_metadata))
  flush.console()
}

write_csv(all_metadata, "raw_sra_metadata.csv", progress = FALSE)
message("Data saved successfully: raw_sra_metadata.csv")
flush.console()
EOF

# Run the R script and ensure console output is displayed immediately
Rscript hmp_16s_query.R | tee real_time_log.txt

echo "Script execution completed."
