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
SRA_FILES="${WORKDIR}/fastq_data/${ACCESSION}*.sra"

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

echo "🔹 Downloading FASTQ files for \$ACCESSION using provider: \$PROVIDER"

# Download FASTQ files
fastq-dl --accession "\$ACCESSION" \
         --provider "\$PROVIDER" \
         --cpus 4 \
         --prefix "\$ACCESSION" \
         --outdir "\${WORKDIR}/fastq_data"

# Check if download was successful
if [[ \$? -ne 0 ]]; then
    echo "❌ ERROR: Failed to download FASTQ files for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Verify that at least one FASTQ exists
shopt -s nullglob
ALL_FASTQS=("\${WORKDIR}/fastq_data"/*.fastq "\${WORKDIR}/fastq_data"/*.fastq.gz)

if [ "\${#ALL_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No FASTQ files found for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Ensure all FASTQs are gzipped
for FILE in "\${WORKDIR}/fastq_data"/*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        if ! gzip "\$FILE"; then
            echo "❌ ERROR: Failed to gzip \$FILE" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
            exit 1
        fi
    fi
done

shopt -s nullglob
ALL_GZ_FASTQS=("\${WORKDIR}/fastq_data"/*.fastq.gz)

if [ "\${#ALL_GZ_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No gzipped FASTQ files found for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
    exit 1
fi

# Move metadata file and append to all_fastq_run_info.tsv
METADATA_FILENAME="\${ACCESSION}-run-info.tsv"
ORIG_METADATA_FILE="\${WORKDIR}/fastq_data/\${METADATA_FILENAME}"
ALL_METADATA_FILE="\${WORKDIR}/metadata/all_fastq_run_info.tsv"

if [[ -f "\$ORIG_METADATA_FILE" ]]; then
    echo "🔹 Moving metadata file to metadata/ folder."
    if ! mv "\$ORIG_METADATA_FILE" "\${WORKDIR}/metadata/"; then
        echo "❌ ERROR: Failed to move metadata file for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi
    
    echo "🔹 Appending metadata to all_fastq_run_info.tsv with locking."
    (
        flock -x 200
        if ! cat "\${WORKDIR}/metadata/\${METADATA_FILENAME}" >> "\$ALL_METADATA_FILE"; then
            echo "❌ ERROR: Failed to append metadata for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
            exit 1
        fi
        rm -f "\${LOCK_FILE}"  # Clean up the lock file
    ) 200>"\${LOCK_FILE}"
    
    echo "🔹 Removing individual metadata file."
    if ! rm -f "\${WORKDIR}/metadata/\${METADATA_FILENAME}"; then
        echo "❌ ERROR: Failed to remove metadata file for \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.err"
        exit 1
    fi
else
    echo "⚠️ WARNING: No metadata file found for \$ACCESSION."
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

echo "✅ Successfully downloaded and gzipped FASTQ files for \$ACCESSION."

# Clean up logs if everything is fine
rm -f "$JOB_SCRIPT"
rm -f "${WORKDIR}/logs/${ACCESSION}.out" "${WORKDIR}/logs/${ACCESSION}.err"
EOF

    chmod +x "$JOB_SCRIPT"
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
    [[ -z "$ACCESSION" ]] && continue

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
