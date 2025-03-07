#!/bin/bash
#SBATCH --job-name=main_DL
#SBATCH --output=%x_logs/slurm_%A.out
#SBATCH --error=%x_logs/slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

#######################################
# Configuration
#######################################
MAX_JOBS=50  # Maximum number of concurrent jobs
WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
LOCK_FILE="${WORKDIR}/checkpoint.lock"
ACCESSIONS_FILE="${WORKDIR}/run_accessions.txt"
COMBINED_METADATA="${WORKDIR}/metadata/combined_metadata.tsv"

# Ensure we are working in the right directory
cd "$WORKDIR" || { echo "❌ ERROR: Unable to change to working directory $WORKDIR"; exit 1; }

# Create necessary directories
mkdir -p "${WORKDIR}/jobs" "${WORKDIR}/logs" "${WORKDIR}/fastq_data" "${WORKDIR}/metadata"

# Ensure checkpoint file exists
touch "${CHECKPOINT_FILE}"

#######################################
# Function to Submit a Job for an Accession
#######################################
submit_job() {
    local ACCESSION=$1
    local JOB_SCRIPT="${WORKDIR}/jobs/download_${ACCESSION}.sh"

    # Trim leading/trailing whitespace from the accession
    ACCESSION=$(echo "$ACCESSION" | xargs)

    # Debug: Log the accession being processed
    echo "DEBUG: Submitting job for accession: $ACCESSION" >> "${WORKDIR}/logs/debug.log"

    # Create the job script
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=${WORKDIR}/logs/${ACCESSION}.out
#SBATCH --error=${WORKDIR}/logs/${ACCESSION}.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# Configuration
WORKDIR="$WORKDIR"
CHECKPOINT_FILE="$CHECKPOINT_FILE"
LOCK_FILE="$LOCK_FILE"
SRA_FILES="${WORKDIR}/fastq_data/${ACCESSION}.sra"
FASTQ_DIR="${WORKDIR}/fastq_data"
METADATA_DIR="${WORKDIR}/metadata"
COMBINED_METADATA="$COMBINED_METADATA"
ACCESSION="$ACCESSION"

# Determine the provider based on the accession prefix
ACCESSION_PREFIX=\${ACCESSION:0:3}
echo "DEBUG: Accession = \$ACCESSION, Prefix = \$ACCESSION_PREFIX" >> "\${WORKDIR}/logs/\${ACCESSION}.out"

if [[ "\$ACCESSION_PREFIX" =~ ^(SRR|SRX|SRS|SRP)$ ]]; then
    PROVIDER="sra"
elif [[ "\$ACCESSION_PREFIX" =~ ^(ERR|ERX|ERS|ERP|DRR|DRX|DRS|DRP)$ ]]; then
    PROVIDER="ena"
else
    echo "❌ ERROR: Unknown accession type: \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

echo "DEBUG: Using provider: \$PROVIDER for accession: \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.out"

if [[ "\$PROVIDER" == "sra" ]]; then
    # Step 1: Prefetch
    echo "🔹 Prefetching SRA file for \$ACCESSION"
    prefetch "\$ACCESSION" --output-file "\$SRA_FILES" || { echo "❌ ERROR: prefetch failed"; exit 1; }
    sleep 3  # Throttle after prefetch

    # Step 2: Convert to FASTQ
    echo "🔹 Converting SRA to FASTQ"
    fasterq-dump "\$ACCESSION" \
        --outdir "\$FASTQ_DIR" \
        --threads 4 \
        --mem 8G \
        --split-3 || { echo "❌ ERROR: fasterq-dump failed"; exit 1; }
    sleep 3  # Throttle after fasterq-dump

    # Step 3: Fetch Metadata with retries
    MAX_RETRIES=3
    RETRY_DELAY=5
    echo "🔹 Fetching metadata (max \$MAX_RETRIES attempts)"
    for ((i=1; i<=\$MAX_RETRIES; i++)); do
        echo "Attempt \$i/3..."
        esearch -db sra -query "\$ACCESSION" | sleep 3 | efetch -format runinfo > "\${METADATA_DIR}/\${ACCESSION}-run-info.csv"
        if [[ \$? -eq 0 && -s "\${METADATA_DIR}/\${ACCESSION}-run-info.csv" ]]; then
            echo "Metadata fetched successfully"
            break
        else
            echo "⚠️ WARNING: Metadata attempt \$i failed"
            sleep \$((RETRY_DELAY * i))
        fi
    done

    # Final check for metadata
    if [[ ! -s "\${METADATA_DIR}/\${ACCESSION}-run-info.csv" ]]; then
        echo "❌ ERROR: All metadata attempts failed"
    fi

elif [[ "\$PROVIDER" == "ena" ]]; then
    # Use fastq-dl for ENA accessions (includes metadata)
    echo "🔹 Downloading FASTQ and metadata for \$ACCESSION"
    fastq-dl --accession "\$ACCESSION" \
             --provider "\$PROVIDER" \
             --cpus 4 \
             --prefix "\$ACCESSION" \
             --outdir "\$FASTQ_DIR" || { echo "❌ ERROR: fastq-dl failed"; exit 1; }

    # Process ENA metadata
    ENA_METADATA="\${FASTQ_DIR}/\${ACCESSION}-run-info.tsv"
    if [[ -f "\$ENA_METADATA" ]]; then
        # Append to combined metadata
        if [[ ! -f "\$COMBINED_METADATA" ]]; then
            head -n 1 "\$ENA_METADATA" > "\$COMBINED_METADATA"
        fi
        tail -n +2 "\$ENA_METADATA" >> "\$COMBINED_METADATA"
        # Move metadata file to metadata directory
        mv "\$ENA_METADATA" "\${METADATA_DIR}/"
        echo "🔹 Metadata processed and moved"
    else
        echo "⚠️ WARNING: ENA metadata file not found"
    fi

else
    echo "❌ ERROR: Unsupported provider: \$PROVIDER" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1  # Important: Exit on unknown provider
fi

# Verify FASTQ files and compress
shopt -s nullglob
for FILE in "\${FASTQ_DIR}"/"\${ACCESSION}"*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE" || { echo "❌ ERROR: Gzip failed"; }
    fi
done

# Check for gzipped FASTQs
ALL_GZ_FASTQS=("\${FASTQ_DIR}"/"\${ACCESSION}"*.fastq.gz)
[[ \${#ALL_GZ_FASTQS[@]} -eq 0 ]] && { echo "❌ ERROR: No FASTQ files"; exit 1; }

# Record success in checkpoint
flock -x "\$LOCK_FILE" -c "echo \$ACCESSION >> \$CHECKPOINT_FILE"

# Cleanup SRA files
if ls "\$SRA_FILES" 1>/dev/null 2>&1; then
    rm -f "\$SRA_FILES" || { echo "❌ ERROR: Failed to delete SRA files"; exit 1; }
fi

# Cleanup job script and logs on success
rm -f "$JOB_SCRIPT" "\${WORKDIR}/logs/\${ACCESSION}.out" "\${WORKDIR}/logs/\${ACCESSION}.err"

echo "✅ Successfully processed \$ACCESSION"
EOF

    chmod +x "$JOB_SCRIPT"
    sbatch "$JOB_SCRIPT" || { echo "❌ ERROR: Failed to submit job"; exit 1; }
    echo "✅ Submitted job for ${ACCESSION}"
}

#######################################
# Main Script Logic
#######################################
TOTAL_JOBS=0

while read -r ACCESSION; do
    [[ -z "$ACCESSION" ]] && continue
    ACCESSION=$(echo "$ACCESSION" | xargs)

    if grep -Fxq "$ACCESSION" "$CHECKPOINT_FILE"; then
        echo "⏩ Skipping $ACCESSION"
        continue
    fi

    # Control job concurrency
    while [[ $(squeue -u $USER -h -j -n "fastq_$ACCESSION" | wc -l) -ge $MAX_JOBS ]]; do
        sleep 60
    done

    submit_job "$ACCESSION"
    TOTAL_JOBS=$((TOTAL_JOBS + 1))
done < "$ACCESSIONS_FILE"

echo "🎉 All $TOTAL_JOBS jobs submitted!"
