#!/bin/bash
#SBATCH --time=00:30:00   # Adjust as needed
#SBATCH --mem=2G          # Adjust memory as needed
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

# Get the study ID from the command line
STUDY_ID=$1
echo "🔍 Processing study: $STUDY_ID"

# Define the output file and a separate lock file
OUTPUT_FILE="run_accessions.txt"
LOCK_FILE="run_accessions.lock"

# Acquire an exclusive lock (via file descriptor 200)
exec 200>"$LOCK_FILE"
flock -x 200   # Blocks until lock is available

# ---[ Start of Critical Section ]---
if [[ $STUDY_ID == PRJNA* ]]; then
    # NCBI SRA query
    esearch -db sra -query "$STUDY_ID" \
      | efetch -format runinfo \
      | cut -d ',' -f 1 \
      | tail -n +2 \
      >> "$OUTPUT_FILE"

elif [[ $STUDY_ID == PRJEB* || $STUDY_ID == ERP* ]]; then
    # ENA query
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY_ID&result=read_run&fields=run_accession&format=tsv" \
      | tail -n +2 \
      >> "$OUTPUT_FILE"

else
    echo "❌ Unknown study type: $STUDY_ID (Skipping)"
fi

echo "✅ Finished processing $STUDY_ID"
