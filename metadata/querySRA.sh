#!/bin/bash
#OBTAINING METADATA FROM SRA

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
                       "tidyr", "readr")

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

# Warn the user if any package failed to load
if (any(!loaded_packages)) {
    warning("The following packages failed to load: ", 
            paste(required_packages[!loaded_packages], collapse = ", "))
    flush.console()
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
      # Fetch metadata in 'native' format
      metadata_text <- entrez_fetch(db = db_name, web_history = search_results\$web_history, 
                                    rettype = "native", retmode = "text", 
                                    retstart = start, retmax = batch_size)
      
      # Convert text into a tibble
      batch_metadata <- read_delim(metadata_text, delim = "\t", col_types = cols(.default = "c"))
      
      message("Successfully fetched batch: ", start)
      flush.console()
      return(batch_metadata)
      
    }, error = function(e) {
      retries <- retries + 1
      message("Retry ", retries, " for batch at ", start, " due to error: ", e\$message)
      flush.console()
      Sys.sleep(delay_time * retries)  # Increase delay with each retry
    })
  }
  
  warning("Failed to fetch batch starting at ", start, " after ", max_retries, " retries.")
  flush.console()
  return(NULL)
}


batch_indices <- seq(0, total_records - 1, by = 10000)  # Larger batch size for sequential processing

all_metadata <- map_dfr(batch_indices, function(start) {
    message("Fetching records ", start, " to ", min(start + 10000, total_records))
    flush.console()
    result <- fetch_batch(start, 10000, search_results, "sra")
    return(result)
})

if (nrow(all_metadata) != total_records) {
  warning("Mismatch in total records. Expected:", total_records, "Got:", nrow(all_metadata))
  flush.console()
} else {
  message("All records fetched successfully. Total rows: ", nrow(all_metadata))
  flush.console()
}

# Ensure all columns are character type
all_metadata <- all_metadata %>%
  mutate(across(everything(), as.character))

# Identify and log list-columns (if any)
problematic_columns <- names(all_metadata)[sapply(all_metadata, is.list)]
if (length(problematic_columns) > 0) {
    message("Warning: The following columns are lists and may cause formatting issues: ", 
            paste(problematic_columns, collapse = ", "))
    
    # Convert list columns to semicolon-separated strings
    all_metadata <- all_metadata %>%
        mutate(across(where(is.list), ~map_chr(.x, ~paste(.x, collapse = "; "))))
}

# Remove duplicate rows
all_metadata <- distinct(all_metadata)

# Truncate overly long fields (>10,000 characters)
truncate_long_fields <- function(x) {
  ifelse(nchar(x) > 10000, substr(x, 1, 10000), x)
}

all_metadata <- all_metadata %>%
  mutate(across(everything(), truncate_long_fields))

# Check for embedded commas in fields
contains_commas <- sapply(all_metadata, function(col) any(grepl(",", col, fixed = TRUE)))
message("Columns containing commas: ", paste(names(all_metadata)[contains_commas], collapse = ", "))

# Write cleaned metadata to CSV
write_csv(all_metadata, "raw_sra_metadata.csv", na = "", progress = FALSE)

message("Data saved successfully: raw_sra_metadata.csv")
flush.console()
EOF

# Run the R script and ensure console output is displayed immediately
Rscript hmp_16s_query.R

echo "Script execution completed."
