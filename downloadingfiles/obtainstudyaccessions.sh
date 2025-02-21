#!/bin/bash

# Input file containing study IDs
STUDY_FILE="studies.txt"

# Create logs directory if it doesn't exist
mkdir -p logs

# Loop through each study and submit a separate SLURM job
while read -r study_id; do
    # Submit the job and pass the study ID as a parameter
    sbatch --job-name=fetch_${study_id} --output=logs/${study_id}.out --error=logs/${study_id}.err fetch_accessions_slurm.sh "$study_id"
    echo "✅ Submitted job for $study_id"
done < "$STUDY_FILE"

