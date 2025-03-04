#!/bin/bash

# Check if R is installed via Conda, install if missing
if ! command -v R &> /dev/null; then
    echo "R is not installed. Installing R using Conda..."
    conda install -y r-base
fi

# Create the R script to fetch metadata
cat <<EOF > hmp_16s_query.R
# Load necessary libraries
install.packages("rentrez")
install.packages("jsonlite")
install.packages("httr")
library(rentrez)
library(jsonlite)
library(httr)

# Define search term for Human Microbiome Project
search_term <- "Human Microbiome Project[Title] AND 16S[Title] AND Illumina[Title]"

# Search SRA for matching studies
search_results <- entrez_search(db="sra", term=search_term, retmax=500)

# Get details for all found IDs
sra_metadata <- entrez_summary(db="sra", id=search_results\$ids)

# Convert to a dataframe
metadata_list <- lapply(sra_metadata, function(x) data.frame(
  Run = x\$accession,
  Study = x\$study,
  Title = x\$title,
  Platform = x\$platform,
  Strategy = x\$librarystrategy,
  Source = x\$librarysource,
  Layout = x\$librarylayout,
  stringsAsFactors = FALSE
))

# Combine all results into a single dataframe
sra_df <- do.call(rbind, metadata_list)

# Save results to CSV
write.csv(sra_df, "HMP_16S_Illumina.csv", row.names=FALSE)

# Display first few rows
print(head(sra_df))

# Install and connect SQLite
install.packages("RSQLite")
library(RSQLite)
db <- dbConnect(SQLite(), "HMP_16S_Illumina.sqlite")
dbWriteTable(db, "sra_metadata", sra_df, overwrite=TRUE)
dbDisconnect(db)
EOF

# Run the R script
Rscript hmp_16s_query.R
