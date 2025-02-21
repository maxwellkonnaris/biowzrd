#!/bin/bash
#SBATCH --time=00:30:00  # Adjust as needed
#SBATCH --mem=2G         # Adjust memory
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

# Get the study ID from the input argument
STUDY_ID=$1
echo "🔍 Processing study: $STUDY_ID"

# Define output file
OUTPUT_FILE="run_accessions.txt"

# Fetch Run Accessions
if [[ $STUDY_ID == PRJNA* ]]; then
    # NCBI SRA query
    esearch -db sra -query "$STUDY_ID" | efetch -format runinfo | cut -d ',' -f 1 | tail -n +2 >> "$OUTPUT_FILE"

elif [[ $STUDY_ID == PRJEB* || $STUDY_ID == ERP* ]]; then
    # ENA query
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY_ID&result=read_run&fields=run_accession&format=tsv" | tail -n +2 >> "$OUTPUT_FILE"

else
    echo "❌ Unknown study type: $STUDY_ID (Skipping)"
fi

echo "✅ Finished processing $STUDY_ID"

