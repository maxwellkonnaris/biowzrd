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
RUN_INFO_LOCK_FILE="${WORKDIR}/run_info.lock"
ALL_RUN_INFO_FILE="${WORKDIR}/metadata/all_fastq_run_info.tsv"

mkdir -p jobs logs fastq_data metadata

# Ensure these files exist
touch "${CHECKPOINT_FILE}"
touch "${ALL_RUN_INFO_FILE}"

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
RUN_INFO_LOCK_FILE="${RUN_INFO_LOCK_FILE}"
ALL_RUN_INFO_FILE="${ALL_RUN_INFO_FILE}"

FASTQ_DIR="${WORKDIR}/fastq_data"
METADATA_DIR="${WORKDIR}/metadata"

mkdir -p "\${FASTQ_DIR}" "\${METADATA_DIR}"

# Detect provider automatically if accession starts with SRR/ERR/DRR
PROVIDER="ena"
if [[ "${ACCESSION}" == SRR* || "${ACCESSION}" == ERR* || "${ACCESSION}" == DRR* ]]; then
    PROVIDER="sra"
fi

echo "🔹 Downloading ${ACCESSION} from \${PROVIDER}"

# We use --prefix so each run-info file is named uniquely: e.g. "SRS123-run-info.tsv"
for attempt in {1..3}; do
    fastq-dl --accession "${ACCESSION}" \\
             --provider "\${PROVIDER}" \\
             --cpus 4 \\
             --prefix "${ACCESSION}" \\
             --outdir "\${FASTQ_DIR}" \\
        && break

    echo "⚠️ Attempt \$attempt failed for ${ACCESSION}, retrying in 60 seconds..."
    sleep 60
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

########################################
# Move & combine run-info file(s)
########################################
# fastq-dl with --prefix => e.g.  SRS123-run-info.tsv
RUN_INFO_FILE="\${FASTQ_DIR}/${ACCESSION}-run-info.tsv"
if [[ -f "\$RUN_INFO_FILE" ]]; then
    # Move it to metadata
    NEW_RUN_INFO="\${METADATA_DIR}/${ACCESSION}-run-info.tsv"
    mv "\$RUN_INFO_FILE" "\$NEW_RUN_INFO"

    # Append to single master file
    (
       flock -x 300
       cat "\$NEW_RUN_INFO" >> "\$ALL_RUN_INFO_FILE"
    ) 300>"\${RUN_INFO_LOCK_FILE}"

    # Remove now that it's merged
    rm -f "\$NEW_RUN_INFO"
else
    echo "⚠️ No run-info TSV found for ${ACCESSION}. Possibly multi-run or no separate file?"
fi

########################################
# Move metadata JSON(s)
########################################
if compgen -G "\${FASTQ_DIR}/*.metadata.json" > /dev/null; then
    mv "\${FASTQ_DIR}"/*.metadata.json "\${METADATA_DIR}/"
fi

########################################
# Record success in checkpoint
########################################
(
    flock -x 200
    echo "${ACCESSION}" >> "\${CHECKPOINT_FILE}"
) 200>"\${LOCK_FILE}"

########################################
# All done – remove logs
########################################
echo "✅ Finished ${ACCESSION}."

rm -f "logs/${ACCESSION}.out" "logs/${ACCESSION}.err"
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
