#!/bin/bash
#SBATCH --job-name=fasterq_array
#SBATCH --array=1-XX
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --account=open

FASTQ_DIR="./fastq_data"
mkdir -p "$FASTQ_DIR"

SRA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_list.txt)
BASENAME=$(basename "$SRA_FILE" .sra)

echo "[Task $SLURM_ARRAY_TASK_ID] Starting on $SRA_FILE"

# Run fasterq-dump
fasterq-dump "$FASTQ_DIR"/"$SRA_FILE" \
    --outdir "$FASTQ_DIR" \
    --threads 4 \
    --mem 8G \
    --split-3

STATUS=$?

if [ $STATUS -eq 0 ]; then
    echo "[Task $SLURM_ARRAY_TASK_ID] Dump complete. Starting compression."

    # Gzip each output file and track success
    GZIP_SUCCESS=true
    for fq in "$FASTQ_DIR/${BASENAME}"*.fastq; do
        gzip -f "$fq" || GZIP_SUCCESS=false
    done

    if $GZIP_SUCCESS; then
        echo "[Task $SLURM_ARRAY_TASK_ID] Compression complete. Deleting $SRA_FILE."
        rm "$FASTQ_DIR"/"$SRA_FILE"
    else
        echo "[Task $SLURM_ARRAY_TASK_ID] Compression failed for some files. Keeping $SRA_FILE."
    fi
else
    echo "[Task $SLURM_ARRAY_TASK_ID] fasterq-dump failed. Skipping compression and deletion."
fi
