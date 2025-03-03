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
LOCK_FILE="checkpoint.lock"

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

# Detect provider: ENA or SRA
PROVIDER="ena"  # Default to ENA
if [[ "${ACCESSION}" == SRR* || "${ACCESSION}" == ERR* || "${ACCESSION}" == DRR* ]]; then
    PROVIDER="sra"
fi

# Download FASTQ files with retries
echo "🔹 Downloading ${ACCESSION} from \${PROVIDER}"
for attempt in {1..3}; do
    fastq-dl -a ${ACCESSION} --provider \${PROVIDER} --cpus 4 -o \${FASTQ_DIR} && break
    echo "⚠️ Attempt \$attempt failed for ${ACCESSION}, retrying in 60 seconds..."
    sleep 60
done

# Ensure FASTQ files exist before proceeding
if ! ls \${FASTQ_DIR}/${ACCESSION}*.fastq* &>/dev/null; then
    echo "❌ ERROR: Download failed for ${ACCESSION}" >> logs/${ACCESSION}.err
    exit 1
fi

# Move metadata JSON file to metadata folder
mv \${FASTQ_DIR}/${ACCESSION}.metadata.json \${METADATA_DIR}/ 2>/dev/null || echo "⚠️ No metadata file found for ${ACCESSION}"

# Ensure all FASTQ files are gzipped
for FILE in \${FASTQ_DIR}/${ACCESSION}*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE"
    fi
done

# Locking mechanism to prevent race condition in checkpointing
(
    flock -x 200
    echo "${ACCESSION}" >> "$CHECKPOINT_FILE"
) 200>"$LOCK_FILE"

echo "✅ Finished ${ACCESSION}, metadata stored in \${METADATA_DIR}/${ACCESSION}.metadata.json"
EOF

    # Make script executable
    chmod +x "$JOB_SCRIPT"

    # Submit the job
    sbatch "$JOB_SCRIPT"
    echo "✅ Submitted job for ${ACCESSION}"

    COUNT=$((COUNT + 1))
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    # Wait if too many jobs are in the queue
    if [[ $COUNT -ge $BATCH_SIZE ]]; then
        echo "⏳ Waiting for batch of $BATCH_SIZE jobs to complete..."
        while [[ $(squeue --format="%j" -u $USER | grep -c "fastq_") -ge $BATCH_SIZE ]]; do
            sleep 60  # Check job queue every 60 seconds
        done
        COUNT=0  # Reset counter after batch completes
    fi

done < run_accessions.txt

echo "🎉 All $TOTAL_JOBS jobs submitted!"
