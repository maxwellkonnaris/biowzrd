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
BATCH_SIZE=50

WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
LOCK_FILE="${WORKDIR}/checkpoint.lock"

mkdir -p jobs logs fastq_data metadata

# Ensure the checkpoint file exists
touch "${CHECKPOINT_FILE}"

COUNT=0
TOTAL_JOBS=0

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

# Detect provider automatically if accession starts with SRR/ERR/DRR
PROVIDER="ena"
if [[ "${ACCESSION}" == SRR* || "${ACCESSION}" == ERR* || "${ACCESSION}" == DRR* ]]; then
    PROVIDER="sra"
fi

# Convert SRS to SRR before downloading (if needed)
RUN_ACCESSIONS="${ACCESSION}"
if [[ "${ACCESSION}" == SRS* || "${ACCESSION}" == ERS* ]]; then
    RUN_ACCESSIONS=\$(curl -ks "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession" | tail -n +2)

    if [[ -z "\$RUN_ACCESSIONS" ]]; then
        echo "❌ ERROR: No runs found for ${ACCESSION} (likely an SRS accession)" >> "logs/${ACCESSION}.err"
        exit 1
    fi
fi

echo "🔹 Downloading FASTQ files for ${ACCESSION} from \${PROVIDER}"

for SRR in \$RUN_ACCESSIONS; do
    fastq-dl --accession "\$SRR" \\
             --provider "\${PROVIDER}" \\
             --cpus 4 \\
             --outdir "\${FASTQ_DIR}"
done

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
