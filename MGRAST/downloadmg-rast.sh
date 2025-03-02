#!/bin/bash
#SBATCH --job-name=mg_download
#SBATCH --output=logs/mg_download_%A_%a.out
#SBATCH --error=logs/mg_download_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low

# Create logs directory if it doesn't exist
mkdir -p logs

# If the SLURM_ARRAY_TASK_ID is not set, we assume this is the wrapper invocation.
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    # Calculate the total number of metagenome IDs (one per line) in mgrastsamples.txt
    NUM_SAMPLES=$(wc -l < mgrastsamples.txt)
    echo "Submitting job array with ${NUM_SAMPLES} tasks (max 50 concurrently)."
    # Submit this script as a job array with at most 50 concurrent tasks.
    sbatch --array=1-${NUM_SAMPLES}%50 "$0"
    exit 0
fi

# At this point, SLURM_ARRAY_TASK_ID is set, so we are running as an array job.
SAMPLES_FILE="mgrastsamples.txt"

# Extract the metagenome ID corresponding to this task's array index.
mgm_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

if [ -z "$mgm_id" ]; then
    echo "Error: No metagenome ID found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Downloading metagenome ID: ${mgm_id}"

# Download the FASTA file and compress it.
curl -k -L "https://api.mg-rast.org/download/${mgm_id}?file=299.1" | gzip > "${mgm_id}.fasta.gz"

if [ $? -eq 0 ]; then
    echo "Download and compression complete for ${mgm_id}"
else
    echo "Download failed for ${mgm_id}"
    exit 1
fi

