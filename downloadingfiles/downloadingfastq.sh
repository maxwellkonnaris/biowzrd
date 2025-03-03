#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

#######################################
# Configuration
#######################################
BATCH_SIZE=75

WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
LOCK_FILE="${WORKDIR}/checkpoint.lock"

mkdir -p jobs logs fastq_data metadata

# Ensure checkpoint file exists
touch "${CHECKPOINT_FILE}"

COUNT=0
TOTAL_JOBS=0

#######################################
# Function to Determine the Correct Provider
#######################################
get_provider() {
    ACCESSION=$1
    PREFIX=${ACCESSION:0:3}

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
            echo "❌ ERROR: BioProject accessions (PRJ...) must be resolved to Studies (SRP/ERP/DRP) first." >> "logs/$ACCESSION.err"
            exit 1
            ;;
        SAM) 
            echo "❌ ERROR: BioSample accessions (SAM...) must be resolved to Samples (SRS/ERS/DRS) first." >> "logs/$ACCESSION.err"
            exit 1
            ;;
        *)
            echo "❌ ERROR: Unknown accession type: $ACCESSION" >> "logs/$ACCESSION.err"
            exit 1
            ;;
    esac
}

#######################################
# Iterate over run_accessions.txt
#######################################
while read -r ACCESSION; do
    # Skip empty lines
    [[ -z "$ACCESSION" ]] && continue

    # Check if ACCESSION is already done
    if grep -Fxq "${ACCESSION}" "${CHECKPOINT_FILE}"; then
        echo "⏩ Skipping ${ACCESSION}, already processed."
        continue
    fi

    #######################################
    # Create an individual job script
    #######################################
    JOB_SCRIPT="jobs/download_${ACCESSION}.sh"
    PROVIDER=$(get_provider "${ACCESSION}")
    if [[ -z "$PROVIDER" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ⏩ Skipping ${ACCESSION}: Invalid or unsupported accession type." >> logs/skipped_accessions.log
        continue  # Skip processing if accession is invalid
    fi
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=logs/${ACCESSION}.out
#SBATCH --error=logs/${ACCESSION}.err
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# Absolute paths
CHECKPOINT_FILE="${CHECKPOINT_FILE}"
LOCK_FILE="${LOCK_FILE}"
FASTQ_DIR="${WORKDIR}/fastq_data"

mkdir -p "\${FASTQ_DIR}"

echo "🔹 Downloading FASTQ files for ${ACCESSION} using provider: ${PROVIDER}"

fastq-dl --accession "${ACCESSION}" \\
         --provider "${PROVIDER}" \\
         --cpus 4 \\
         --prefix "${ACCESSION}" \\
         --outdir "\${FASTQ_DIR}"

########################################
# Move metadata file and append to all_fastq_run_info.tsv
########################################
METADATA_FILE="${ACCESSION}-run-info.tsv"
ALL_METADATA_FILE="metadata/all_fastq_run_info.tsv"

if [[ -f "\${METADATA_FILE}" ]]; then
    echo "🔹 Moving metadata file to metadata directory."
    mv "\${METADATA_FILE}" metadata/
    
    echo "🔹 Appending metadata to all_fastq_run_info.tsv with locking."
    (
        flock -x 200  # Acquire exclusive lock
        cat "metadata/\${METADATA_FILE}" >> "\${ALL_METADATA_FILE}"
    ) 200>"\${LOCK_FILE}"  # Lock on a dedicated file
    
    echo "🔹 Removing individual metadata file."
    rm -f "metadata/\${METADATA_FILE}"
else
    echo "⚠️ WARNING: No metadata file found for ${ACCESSION}."
fi

########################################
# Verify that at least one FASTQ exists
########################################
shopt -s nullglob
ALL_FASTQS=( "\${FASTQ_DIR}"/*.fastq "\${FASTQ_DIR}"/*.fastq.gz )

if [ "\${#ALL_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No FASTQ files found for ${ACCESSION}" >> "logs/${ACCESSION}.err"
    exit 1
fi

########################################
# Ensure all FASTQs are gzipped
########################################
for FILE in "\${FASTQ_DIR}"/*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE"
    fi
done

# Re-check if gzipped FASTQs exist
shopt -s nullglob
ALL_GZ_FASTQS=( "\${FASTQ_DIR}"/*.fastq.gz )

if [ "\${#ALL_GZ_FASTQS[@]}" -eq 0 ]; then
    echo "❌ ERROR: No gzipped FASTQ files found for ${ACCESSION}, skipping completion record." >> "logs/${ACCESSION}.err"
    exit 1
fi

########################################
# Record success in checkpoint
########################################
(
    flock -x 200
    echo "${ACCESSION}" >> "\${CHECKPOINT_FILE}"
) 200>"\${LOCK_FILE}"

########################################
# Clean up logs if everything is fine
########################################
echo "✅ Successfully downloaded and gzipped FASTQ files for ${ACCESSION}."
rm -f "logs/${ACCESSION}.out" "logs/${ACCESSION}.err"
rm -f "jobs/download_${ACCESSION}.sh"
EOF

    chmod +x "$JOB_SCRIPT"

    # Submit job
    sbatch "$JOB_SCRIPT"
    echo "✅ Submitted job for ${ACCESSION}"

    COUNT=$((COUNT + 1))
    TOTAL_JOBS=$((TOTAL_JOBS + 1))

    #######################################
    # Batch queue limit
    #######################################
    if [[ $COUNT -ge $BATCH_SIZE ]]; then
        echo "⏳ Waiting for batch of $BATCH_SIZE jobs to complete..."
        while [ "$(squeue --format="%j" -u $USER | grep -c "fastq_")" -ge "$BATCH_SIZE" ]; do
            sleep 60
        done
        COUNT=0
    fi

done < run_accessions.txt

echo "🎉 All $TOTAL_JOBS jobs submitted!"
