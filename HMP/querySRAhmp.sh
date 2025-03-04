#!/bin/bash

# Check if R is installed, install if missing
if ! command -v R &> /dev/null; then
    echo "R is not installed. Installing R..."
    conda install -y r-base
fi

# Create an R script
cat <<EOF > hmp_16s_query.R
# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SRAdb", ask=FALSE)
install.packages("RSQLite", dependencies=TRUE)

# Load libraries
library(SRAdb)
library(RSQLite)

# Download the SRA metadata database
sqlfile <- getSRAdbFile()

# Connect to the database
sra_con <- dbConnect(RSQLite::SQLite(), sqlfile)

# Define query for HMP 16S sequencing on Illumina
query <- "
SELECT run.run_accession, run.*, experiment.*, study.*
FROM run
JOIN experiment ON run.experiment_accession = experiment.experiment_accession
JOIN study ON experiment.study_accession = study.study_accession
WHERE study.study_title LIKE 'Human Microbiome Project%'
AND experiment.library_strategy = 'AMPLICON'
AND experiment.platform LIKE '%Illumina%'"

# Execute query
hmp_16s_illumina <- dbGetQuery(sra_con, query)

# Save results to a CSV file
write.csv(hmp_16s_illumina, 'HMP_16S_Illumina.csv', row.names=FALSE)

# Show first few rows
print(head(hmp_16s_illumina))

# Disconnect from database
dbDisconnect(sra_con)
EOF

# Run the R script
Rscript hmp_16s_query.R
