#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Number of jobs to submit at a time
BATCH_SIZE=50

# Use absolute paths where possible
WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
LOCK_FILE="${WORKDIR}/checkpoint.lock"
RUN_INFO_LOCK_FILE="${WORKDIR}/run_info.lock"
ALL_RUN_INFO_FILE="${WORKDIR}/metadata/all_fastq_run_info.tsv"

# Create required directories
mkdir -p jobs logs fastq_data metadata

# Create checkpoint file if it doesn't exist
touch "${CHECKPOINT_FILE}"

# Create the combined run-info file if it doesn't exist
touch "${ALL_RUN_INFO_FILE}"

# Counter for batch submission
COUNT=0
TOTAL_JOBS=0

# Loop through each accession in run_accessions.txt
while read -r ACCESSION; do

    # Skip empty lines
    [[ -z "$ACCESSION" ]] && continue

    # Skip accessions that have already been processed
    if grep -Fxq "${ACCESSION}" "${CHECKPOINT_FILE}"; then
        echo "⏩ Skipping ${ACCESSION}, already processed."
        continue
    fi

    JOB_SCRIPT="jobs/download_${ACCESSION}.sh"

    # Create the individual SLURM job script
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=logs/${ACCESSION}.out
#SBATCH --error=logs/${ACCESSION}.err
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

########################
# Absolute paths
########################
CHECKPOINT_FILE="${CHECKPOINT_FILE}"
LOCK_FILE="${LOCK_FILE}"
RUN_INFO_LOCK_FILE="${RUN_INFO_LOCK_FILE}"
ALL_RUN_INFO_FILE="${ALL_RUN_INFO_FILE}"

FASTQ_DIR="${WORKDIR}/fastq_data"
METADATA_DIR="${WORKDIR}/metadata"

mkdir -p "\${FASTQ_DIR}" "\${METADATA_DIR}"

# Detect provider: ENA or SRA
PROVIDER="ena"  # Default to ENA
if [[ "${ACCESSION}" == SRR* || "${ACCESSION}" == ERR* || "${ACCESSION}" == DRR* ]]; then
    PROVIDER="sra"
fi

echo "🔹 Downloading ${ACCESSION} from \${PROVIDER}"
for attempt in {1..3}; do
    fastq-dl -a ${ACCESSION} --provider "\${PROVIDER}" --cpus 4 -o "\${FASTQ_DIR}" && break
    echo "⚠️ Attempt \$attempt failed for ${ACCESSION}, retrying in 60 seconds..."
    sleep 60
done

# Verify FASTQ files were produced
if ! ls "\${FASTQ_DIR}/${ACCESSION}"*.fastq* &>/dev/null; then
    echo "❌ ERROR: Download failed for ${ACCESSION}" >> "logs/${ACCESSION}.err"
    exit 1
fi

# Gzip any uncompressed FASTQs
for FILE in "\${FASTQ_DIR}/${ACCESSION}"*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE"
    fi
done

########################
# Move & rename run-info
########################
# By default, fastq-dl writes "fastq-run-info.tsv" in the output directory.
# So rename it to avoid overwriting across different accessions.
if [[ -f "\${FASTQ_DIR}/fastq-run-info.tsv" ]]; then
    RUN_INFO_FILE="\${METADATA_DIR}/${ACCESSION}.fastq-run-info.tsv"
    mv "\${FASTQ_DIR}/fastq-run-info.tsv" "\$RUN_INFO_FILE"

    # Append to the combined file (all_fastq_run_info.tsv), with a lock
    (
       flock -x 300
       cat "\$RUN_INFO_FILE" >> "\$ALL_RUN_INFO_FILE"
    ) 300>"\${RUN_INFO_LOCK_FILE}"
fi

########################
# Move metadata JSON
########################
# Move the metadata JSON file to the metadata folder (may not exist for every provider)
if [[ -f "\${FASTQ_DIR}/${ACCESSION}.metadata.json" ]]; then
    mv "\${FASTQ_DIR}/${ACCESSION}.metadata.json" "\${METADATA_DIR}/"
else
    echo "⚠️ No metadata JSON found for ${ACCESSION}"
fi

########################
# Record in checkpoint
########################
(
    flock -x 200
    echo "${ACCESSION}" >> "\${CHECKPOINT_FILE}"
) 200>"\${LOCK_FILE}"

echo "✅ Finished ${ACCESSION}, run-info in \${METADATA_DIR}/${ACCESSION}.fastq-run-info.tsv"
EOF

    chmod +x "$JOB_SCRIPT"

    # Submit the job
    sbatch "$JOB_SCRIPT"
    echo "✅ Submitted job for ${ACCESSION}"

    COUNT=$((COUNT + 1))
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    # If we’ve hit the batch size, wait until the queue has fewer jobs
    if [[ $COUNT -ge $BATCH_SIZE ]]; then
        echo "⏳ Waiting for batch of $BATCH_SIZE jobs to complete..."
        while [[ \$(squeue --format="%j" -u $USER | grep -c "fastq_") -ge $BATCH_SIZE ]]; do
            sleep 60
        done
        COUNT=0
    fi

done < run_accessions.txt

echo "🎉 All $TOTAL_JOBS jobs submitted!"
