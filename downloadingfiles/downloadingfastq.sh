#!/bin/bash
#SBATCH --job-name=fastq_download
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
#SBATCH --cpus-per-task=4  # Increase threads for fasterq-dump
#SBATCH --ntasks=1

# Configuration
WORKDIR="$WORKDIR"
CHECKPOINT_FILE="$CHECKPOINT_FILE"
LOCK_FILE="$LOCK_FILE"
SRA_FILES="${WORKDIR}/fastq_data/${ACCESSION}.sra"
FASTQ_DIR="${WORKDIR}/fastq_data"
METADATA_DIR="${WORKDIR}/metadata"

# Pass the accession as a variable
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

# Download and process based on provider
if [[ "\$PROVIDER" == "sra" ]]; then
    # Step 1: Prefetch the SRA file
    echo "🔹 Prefetching SRA file for \$ACCESSION"
    prefetch "\$ACCESSION" --output-file "\$SRA_FILES"

    if [[ \$? -ne 0 ]]; then
        echo "❌ ERROR: prefetch failed for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi

    # Step 2: Convert SRA to FASTQ using fasterq-dump
    echo "🔹 Converting SRA to FASTQ for \$ACCESSION"
    fasterq-dump "\$ACCESSION" \
                 --outdir "\$FASTQ_DIR" \
                 --threads 4 \
                 --mem 8G \
                 --split-3

    if [[ \$? -ne 0 ]]; then
        echo "❌ ERROR: fasterq-dump failed for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi

    # Step 3: Fetch metadata using fastq-dl (metadata only)
    echo "🔹 Fetching metadata for \$ACCESSION using fastq-dl"
    fastq-dl --accession "\$ACCESSION" \
             --provider "\$PROVIDER" \
             --cpus 4 \
             --prefix "\$ACCESSION" \
             --outdir "\$METADATA_DIR" \
             --metadata-only

    if [[ \$? -ne 0 ]]; then
        echo "⚠️ WARNING: Failed to fetch metadata for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    else
        echo "🔹 Metadata saved to \${METADATA_DIR}/\${ACCESSION}-run-info.tsv"
    fi

elif [[ "\$PROVIDER" == "ena" ]]; then
    # Use fastq-dl for ENA accessions (includes metadata by default)
    echo "🔹 Downloading FASTQ files for \$ACCESSION using provider: \$PROVIDER"
    fastq-dl --accession "\$ACCESSION" \
             --provider "\$PROVIDER" \
             --cpus 4 \
             --prefix "\$ACCESSION" \
             --outdir "\$FASTQ_DIR"

    if [[ \$? -ne 0 ]]; then
        echo "❌ ERROR: fastq-dl failed for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi
else
    echo "❌ ERROR: Unsupported provider: \$PROVIDER" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Verify that at least one FASTQ exists
shopt -s nullglob
ALL_FASTQS=("\${FASTQ_DIR}"/*.fastq "\${FASTQ_DIR}"/*.fastq.gz)

if [ "\${#ALL_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No FASTQ files found for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Ensure all FASTQs are gzipped
for FILE in "\${FASTQ_DIR}"/*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        if ! gzip "\$FILE"; then
            echo "❌ ERROR: Failed to gzip \$FILE" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
            exit 1
        fi
    fi
done

shopt -s nullglob
ALL_GZ_FASTQS=("\${FASTQ_DIR}"/*.fastq.gz)

if [ "\${#ALL_GZ_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No gzipped FASTQ files found for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Record success in checkpoint
(
    flock -x 200
    echo "\$ACCESSION" >> "\$CHECKPOINT_FILE"
    rm -f "\${LOCK_FILE}"  # Clean up the lock file
) 200>"\${LOCK_FILE}"

# Delete all .sra files associated with the accession after FASTQ files are successfully downloaded
if ls \$SRA_FILES 1> /dev/null 2>&1; then
    echo "🔹 Removing leftover .sra files: \$SRA_FILES"
    if ! rm -f \$SRA_FILES; then
        echo "❌ ERROR: Failed to delete .sra files for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi
else
    echo "🔹 No .sra files found to delete."
fi

echo "✅ Successfully downloaded and processed FASTQ files for \$ACCESSION."

# Clean up logs if everything is fine
rm -f "$JOB_SCRIPT"
rm -f "${WORKDIR}/logs/${ACCESSION}.out" "${WORKDIR}/logs/${ACCESSION}.err"
EOF

    # Make the job script executable
    chmod +x "$JOB_SCRIPT"

    # Submit the job
    if ! sbatch "$JOB_SCRIPT"; then
        echo "❌ ERROR: Failed to submit job for ${ACCESSION}" >> "${WORKDIR}/logs/${ACCESSION}.err"
        exit 1
    fi
    echo "✅ Submitted job for ${ACCESSION}"
}

#######################################
# Main Script Logic
#######################################
TOTAL_JOBS=0

while read -r ACCESSION; do
    # Skip empty lines
    [[ -z "$ACCESSION" ]] && { echo "⚠️ WARNING: Skipping empty accession."; continue; }

    # Trim leading/trailing whitespace from the accession
    ACCESSION=$(echo "$ACCESSION" | xargs)

    # Log the accession being processed
    echo "Processing accession: $ACCESSION" >> "${WORKDIR}/logs/debug.log"

    # Check if ACCESSION is already done
    if grep -Fxq "${ACCESSION}" "${CHECKPOINT_FILE}"; then
        echo "⏩ Skipping ${ACCESSION}, already processed."
        continue
    fi

    # Wait until the number of running jobs is below MAX_JOBS
    while [ "$(squeue -u $USER --format="%j" | grep -c "fastq_")" -ge "$MAX_JOBS" ]; do
        sleep 60  # Wait for 1 minute before checking again
    done

    # Submit a job for the accession
    submit_job "$ACCESSION"
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

done < "$ACCESSIONS_FILE"

echo "🎉 All $TOTAL_JOBS jobs submitted!"
