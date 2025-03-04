#!/bin/bash

# Check if R is installed via Conda, install if missing
if ! command -v R &> /dev/null; then
    echo "R is not installed. Installing R using Conda..."
    conda install -y r-base
fi

# Define the SRA database file path
SRA_DB="SRAmetadb.sqlite"

# Check if the SRA database already exists
if [ ! -f "$SRA_DB" ]; then
    echo "SRAmetadb.sqlite not found. Downloading..."
    wget https://s3.amazonaws.com/starbuck1/sradb/SRAmetadb.sqlite -O SRAmetadb.sqlite
fi

# Create an R script
cat <<EOF > hmp_16s_query.R
# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SRAdb", ask=FALSE)
install.packages("RSQLite", dependencies=TRUE)

# Load libraries
library(SRAdb)
library(RSQLite)

# Define the database path
sqlfile <- "SRAmetadb.sqlite"

# Connect to the database
sra_con <- dbConnect(RSQLite::SQLite(), sqlfile)

# Use getSRA() to retrieve metadata for HMP 16S Illumina sequencing
hmp_data <- getSRA(search_terms = "Human Microbiome Project", con = sra_con)

# Filter for only 16S sequencing and Illumina platform
hmp_16s_illumina <- subset(hmp_data, library_strategy == "AMPLICON" & grepl("Illumina", platform))

# Save results to a CSV file
write.csv(hmp_16s_illumina, 'HMP_16S_Illumina.csv', row.names=FALSE)

# Show first few rows
print(head(hmp_16s_illumina))

# Disconnect from database
dbDisconnect(sra_con)
EOF

# Run the R script
Rscript hmp_16s_query.R
