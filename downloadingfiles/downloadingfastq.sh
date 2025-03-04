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
BATCH_SIZE=50  # Number of accessions per job
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
# Function to Determine the Correct Provider
#######################################
get_provider() {
    local ACCESSION=$1
    local PREFIX=${ACCESSION:0:3}

    case "$PREFIX" in
        SRR) echo "sra" ;;  # NCBI SRA Run
        ERR|DRR) echo "ena" ;;  # ENA/DRR Run (European/Japanese)
        SRX) echo "sra" ;;  # SRA Experiment
        ERX|DRX) echo "ena" ;;  # ENA/DRX Experiment
        SRS) echo "sra" ;;  # SRA Sample
        ERS|DRS) echo "ena" ;;  # ENA Sample
        SRP) echo "sra" ;;  # SRA Study
        ERP|DRP) echo "ena" ;;  # ENA/DRP Study
        PRJ)
            echo "❌ ERROR: BioProject accessions (PRJ...) must be resolved to Studies (SRP/ERP/DRP) first." >> "${WORKDIR}/logs/${ACCESSION}.err"
            exit 1
            ;;
        SAM)
            echo "❌ ERROR: BioSample accessions (SAM...) must be resolved to Samples (SRS/ERS/DRS) first." >> "${WORKDIR}/logs/${ACCESSION}.err"
            exit 1
            ;;
        *)
            echo "❌ ERROR: Unknown accession type: $ACCESSION" >> "${WORKDIR}/logs/${ACCESSION}.err"
            exit 1
            ;;
    esac
}

#######################################
# Function to Process a Batch of Accessions
#######################################
process_batch() {
    local BATCH=("$@")
    local JOB_SCRIPT="${WORKDIR}/jobs/download_batch_$(date +%s).sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_batch
#SBATCH --output=${WORKDIR}/logs/batch_%j.out
#SBATCH --error=${WORKDIR}/logs/batch_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

# Set working directory
cd "${WORKDIR}" || { echo "❌ ERROR: Unable to change to working directory ${WORKDIR}"; exit 1; }

FASTQ_DIR="${WORKDIR}/fastq_data"

# Function to download a single accession
download_accession() {
    local ACCESSION=\$1
    local PROVIDER=\$(get_provider "\$ACCESSION")

    echo "🔹 Downloading FASTQ files for \${ACCESSION} using provider: \${PROVIDER}"

    # Download FASTQ files
    fastq-dl --accession "\${ACCESSION}" \\
             --provider "\${PROVIDER}" \\
             --cpus 2 \\
             --prefix "\${ACCESSION}" \\
             --outdir "\${FASTQ_DIR}"

    # Check if download was successful
    if [[ \$? -ne 0 ]]; then
        echo "❌ ERROR: Failed to download FASTQ files for \${ACCESSION}" >> "${WORKDIR}/logs/\${ACCESSION}.err"
        return 1
    fi

    # Move metadata file and append to all_fastq_run_info.tsv
    METADATA_FILENAME="\${ACCESSION}-run-info.tsv"
    ORIG_METADATA_FILE="\${FASTQ_DIR}/\${METADATA_FILENAME}"
    ALL_METADATA_FILE="${WORKDIR}/metadata/all_fastq_run_info.tsv"

    if [[ -f "\${ORIG_METADATA_FILE}" ]]; then
        echo "🔹 Moving metadata file to metadata/ folder."
        mv "\${ORIG_METADATA_FILE}" "${WORKDIR}/metadata/"
        
        echo "🔹 Appending metadata to all_fastq_run_info.tsv with locking."
        (
            flock -x 200  # Acquire exclusive lock
            cat "${WORKDIR}/metadata/\${METADATA_FILENAME}" >> "\${ALL_METADATA_FILE}"
        ) 200>"${LOCK_FILE}"
        
        echo "🔹 Removing individual metadata file."
        rm -f "${WORKDIR}/metadata/\${METADATA_FILENAME}"
    else
        echo "⚠️ WARNING: No metadata file found for \${ACCESSION}."
    fi

    # Verify that at least one FASTQ exists
    shopt -s nullglob
    ALL_FASTQS=( "\${FASTQ_DIR}"/*.fastq "\${FASTQ_DIR}"/*.fastq.gz )

    if [ "\${#ALL_FASTQS[@]}" -eq 0 ]; then
        echo "❌ ERROR: No FASTQ files found for \${ACCESSION}" >> "${WORKDIR}/logs/\${ACCESSION}.err"
        return 1
    fi

    # Ensure all FASTQs are gzipped
    for FILE in "\${FASTQ_DIR}"/*.fastq; do
        if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
            echo "🔹 Gzipping: \$FILE"
            gzip "\$FILE"
        fi
    done

    shopt -s nullglob
    ALL_GZ_FASTQS=( "\${FASTQ_DIR}"/*.fastq.gz )

    if [ "\${#ALL_GZ_FASTQS[@]}" -eq 0 ]; then
        echo "❌ ERROR: No gzipped FASTQ files found for \${ACCESSION}, skipping completion record." >> "${WORKDIR}/logs/\${ACCESSION}.err"
        return 1
    fi

    # Record success in checkpoint
    (
        flock -x 200
        echo "\${ACCESSION}" >> "${CHECKPOINT_FILE}"
    ) 200>"${LOCK_FILE}"

    # Delete .sra if present
    SRA_FILE="\${FASTQ_DIR}/\${ACCESSION}.sra"
    if [[ -f "\${SRA_FILE}" ]]; then
        echo "🔹 Removing leftover .sra file: \${SRA_FILE}"
        rm -f "\${SRA_FILE}"
    fi

    echo "✅ Successfully downloaded and gzipped FASTQ files for \${ACCESSION}."
}

# Export functions and variables for GNU Parallel
export -f download_accession get_provider
export WORKDIR FASTQ_DIR LOCK_FILE CHECKPOINT_FILE

# Process accessions in parallel
echo "🔹 Processing batch of ${#BATCH[@]} accessions in parallel..."
parallel -j 4 download_accession ::: "${BATCH[@]}"

# Clean up logs if everything is fine
rm -f "${JOB_SCRIPT}"
EOF

    chmod +x "$JOB_SCRIPT"
    sbatch "$JOB_SCRIPT"
    echo "✅ Submitted batch job for ${#BATCH[@]} accessions."
}

#######################################
# Main Script Logic
#######################################
COUNT=0
TOTAL_JOBS=0
BATCH=()

while read -r ACCESSION; do
    # Skip empty lines
    [[ -z "$ACCESSION" ]] && continue

    # Check if ACCESSION is already done
    if grep -Fxq "${ACCESSION}" "${CHECKPOINT_FILE}"; then
        echo "⏩ Skipping ${ACCESSION}, already processed."
        continue
    fi

    # Add accession to batch
    BATCH+=("$ACCESSION")
    COUNT=$((COUNT + 1))
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    # Submit batch if size is reached
    if [[ $COUNT -ge $BATCH_SIZE ]]; then
        process_batch "${BATCH[@]}"
        BATCH=()
        COUNT=0
    fi

done < "$ACCESSIONS_FILE"

# Submit remaining accessions in the last batch
if [[ ${#BATCH[@]} -gt 0 ]]; then
    process_batch "${BATCH[@]}"
fi

echo "🎉 All $TOTAL_JOBS accessions submitted in batches!"
