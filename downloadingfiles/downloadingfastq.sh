#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Number of jobs to submit at a time
BATCH_SIZE=50
CHECKPOINT_FILE="completed_accessions.txt"

# Create required directories
mkdir -p jobs logs fastq_data metadata

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"

# Counter for batch submission
COUNT=0
TOTAL_JOBS=0

# Loop through each accession in run_accessions.txt
while read -r ACCESSION; do

    # Skip accessions that have already been processed
    if grep -Fxq "${ACCESSION}" "$CHECKPOINT_FILE"; then
        echo "⏩ Skipping ${ACCESSION}, already processed."
        continue
    fi

    JOB_SCRIPT="jobs/download_${ACCESSION}.sh"

    # Create individual SLURM job script for each accession
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=logs/${ACCESSION}.out
#SBATCH --error=logs/${ACCESSION}.err
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

FASTQ_DIR="fastq_data"
METADATA_DIR="metadata/"

mkdir -p \${FASTQ_DIR} \${METADATA_DIR}

# Check provider type: ENA or SRA
PROVIDER="ena"  # Default to ENA

if [[ "${ACCESSION}" == SRR* || "${ACCESSION}" == ERR* || "${ACCESSION}" == DRR* ]]; then
    PROVIDER="sra"
fi

# Download FASTQ files with 4 threads and retries
fastq-dl -a ${ACCESSION} --provider \${PROVIDER} --cpus 4 -o \${FASTQ_DIR}

# Move metadata JSON file to metadata folder
mv \${FASTQ_DIR}/${ACCESSION}.metadata.json \${METADATA_DIR}/

# Ensure all FASTQ files are gzipped
for FILE in \${FASTQ_DIR}/${ACCESSION}*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE"
    fi
done

# Mark this accession as completed
echo "${ACCESSION}" >> "$CHECKPOINT_FILE"

echo "✅ Finished ${ACCESSION}, metadata stored in \${METADATA_DIR}/${ACCESSION}.metadata.json"
EOF

    # Make script executable
    chmod +x "$JOB_SCRIPT"

    # Submit the job
    sbatch "$JOB_SCRIPT"
    echo "✅ Submitted job for ${ACCESSION}"

    COUNT=$((COUNT + 1))
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    # If 50 jobs have been submitted, wait for them to finish before submitting more
    if [[ $COUNT -ge $BATCH_SIZE ]]; then
        echo "⏳ Waiting for batch of $BATCH_SIZE jobs to complete..."
        while [[ $(squeue -u $USER | grep -c "fastq_") -ge $BATCH_SIZE ]]; do
            sleep 60  # Check job queue every 60 seconds
        done
        COUNT=0  # Reset counter after batch completes
    fi

done < run_accessions.txt

echo "🎉 All $TOTAL_JOBS jobs submitted!"

