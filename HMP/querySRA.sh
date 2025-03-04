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

# Create the R script dynamically with the user-provided query
cat <<EOF > hmp_16s_query.R
# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load required libraries, install if missing
required_packages <- c("rentrez", "dplyr", "stringr", "purrr", "xml2", "tidyr", "readr", "foreach", "doParallel")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Parallel Processing Setup
num_cores <- parallel::detectCores() - 1  # Use all available cores minus one
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Define search parameters (Injected from Bash script)
search_query <- "$SEARCH_QUERY"
db_name <- "sra"
batch_size <- 500  # Adjustable
delay_time <- 1  # Time delay to avoid API limits
max_retries <- 3  # Maximum number of retries for failed batches

# Print the query for debugging
cat("Search Query:", search_query, "\n")

# Perform initial search
search_results <- tryCatch({
  result <- entrez_search(db = db_name, term = search_query, retmax = 0, use_history = TRUE)
  cat("API Response Received.\n")
  result
}, error = function(e) {
  stop("NCBI API error: ", e$message)
})

# Check if the search returned results
if (is.null(search_results$count)) {
  cat("No records found. Check your query or API key.\n")
} else {
  # Convert total_records to a string to avoid errors in cat()
  total_records <- as.character(search_results$count)
  cat("Total Records Found:", total_records, "\n")
}

# If no records are found, test with a simpler query
if (exists("total_records") && total_records == "0") {
  cat("Testing with a simpler query: 'Human Microbiome Project'\n")
  test_query <- "Human Microbiome Project"
  test_results <- tryCatch({
    result <- entrez_search(db = db_name, term = test_query, retmax = 0, use_history = TRUE)
    cat("API Response for Test Query Received.\n")
    result
  }, error = function(e) {
    stop("NCBI API error: ", e$message)
  })
  
  if (is.null(test_results$count)) {
    cat("No records found even with a simpler query. Check your API key or network connection.\n")
  } else {
    cat("Test Query Records Found:", as.character(test_results$count), "\n")
  }
}

# Function to fetch and parse a batch of records
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

# Parallel Batch Processing
batch_indices <- seq(0, total_records, by = batch_size)

all_metadata <- foreach(start = batch_indices, .combine = bind_rows, .packages = c("rentrez", "dplyr", "xml2", "tidyr")) %dopar% {
  cat("Fetching records", start, "to", min(start + batch_size, total_records), "\\n")
  fetch_batch(start, batch_size, search_results, db_name)
}

# Stop parallel processing
parallel::stopCluster(cl)

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
