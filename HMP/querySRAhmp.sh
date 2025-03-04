#!/bin/bash

# Check if R is installed via Conda, install if missing
if ! command -v R &> /dev/null; then
    echo "R is not installed. Installing R using Conda..."
    conda install -y r-base
fi

# Create an R script to download SRAmetadb.sqlite and query it
cat <<EOF > hmp_16s_query.R
# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SRAdb", ask=FALSE)
install.packages("RSQLite", dependencies=TRUE)

# Load libraries
library(SRAdb)
library(RSQLite)

# Define the database path
sqlfile <- "SRAmetadb.sqlite"

# Download the SRA database if it doesn't exist
if (!file.exists(sqlfile)) {
    message("Downloading SRAmetadb.sqlite, this may take some time...")
    sqlfile <- getSRAdbFile()
}

# Verify that the file was downloaded
if (!file.exists(sqlfile)) {
    stop("ERROR: Failed to download SRAmetadb.sqlite. Please check your network connection or try again later.")
}

# Connect to the SRA database
sra_con <- dbConnect(RSQLite::SQLite(), sqlfile)

# Check available tables
tables <- dbListTables(sra_con)
print(tables)  # Print available tables for debugging

# Determine the correct table name for runs
if ("run" %in% tables) {
    table_name <- "run"
} else if ("sra_run" %in% tables) {
    table_name <- "sra_run"
} else {
    stop("Neither 'run' nor 'sra_run' table found in the database.")
}

# Construct SQL query dynamically
query <- sprintf("
SELECT %s.run_accession, %s.*, experiment.*, study.*
FROM %s
JOIN experiment ON %s.experiment_accession = experiment.experiment_accession
JOIN study ON experiment.study_accession = study.study_accession
WHERE study.study_title LIKE 'Human Microbiome Project%%'
AND experiment.library_strategy = 'AMPLICON'
AND experiment.platform LIKE '%%Illumina%%'", table_name, table_name, table_name, table_name)

# Execute the query
hmp_16s_illumina <- dbGetQuery(sra_con, query)

# Save results to CSV
write.csv(hmp_16s_illumina, 'HMP_16S_Illumina.csv', row.names=FALSE)

# Show first few rows
print(head(hmp_16s_illumina))

# Disconnect from database
dbDisconnect(sra_con)
EOF

# Run the R script
Rscript hmp_16s_query.R
