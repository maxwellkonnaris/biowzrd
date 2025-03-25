#!/bin/bash
#SBATCH --job-name=mg_download_wrapper
#SBATCH --output=logs/mg_download_wrapper.out
#SBATCH --error=logs/mg_download_wrapper.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low

# Create logs directory if it doesn't exist
mkdir -p logs

# Check if this is the wrapper job (SLURM_ARRAY_TASK_ID is not set)
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    # Calculate total number of metagenome IDs in mgrastsamples.txt
    NUM_SAMPLES=$(wc -l < mgrastsamples.txt)
    echo "Submitting job array with ${NUM_SAMPLES} tasks (max 50 concurrently)."
    # Submit this script as a job array with --time=12:00:00 for each download job
    sbatch --time=00:30:00 --array=1-${NUM_SAMPLES}%50 "$0"
    exit 0
fi

# --------- Array Task: Download a specific metagenome ---------

# File containing the list of metagenome IDs
SAMPLES_FILE="mgrastsamples.txt"

# Extract the metagenome ID corresponding to this task's array index.
mgm_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")
if [ -z "$mgm_id" ]; then
    echo "Error: No metagenome ID found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Checkpoint: if the output file already exists and is non-empty, skip this task.
OUTPUT_FILE="${mgm_id}.fasta.gz"
if [ -s "$OUTPUT_FILE" ]; then
    echo "Checkpoint: ${OUTPUT_FILE} already exists. Skipping download for ${mgm_id}."
    exit 0
fi

echo "Downloading metagenome ID: ${mgm_id}"

# Download the FASTA file and compress it on the fly.
curl -k -L "https://api.mg-rast.org/download/${mgm_id}?file=299.1" | gzip > "$OUTPUT_FILE"

if [ $? -eq 0 ]; then
    echo "Download and compression complete for ${mgm_id}"
else
    echo "Download failed for ${mgm_id}"
    exit 1
fi

