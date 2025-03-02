#!/bin/bash
#SBATCH --job-name=mg_download
#SBATCH --output=logs/mg_download_%A_%a.out
#SBATCH --error=logs/mg_download_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low
# Note: The array range (e.g., --array=1-<NUM_SAMPLES>%50) is provided when submitting via the wrapper script.

# Create a logs directory if it doesn't exist.
mkdir -p logs

# File containing the list of metagenome IDs.
SAMPLES_FILE="mgrastsamples.txt"

# Extract the metagenome ID for the current SLURM array task.
mgm_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")

if [ -z "$mgm_id" ]; then
    echo "Error: No metagenome ID found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Downloading metagenome ID: ${mgm_id}"

# Download the FASTA file from MG-RAST and compress it using gzip.
curl -k -L "https://api.mg-rast.org/download/${mgm_id}?file=299.1" | gzip > "${mgm_id}.fasta.gz"

if [ $? -eq 0 ]; then
    echo "Download and compression complete for ${mgm_id}"
else
    echo "Download failed for ${mgm_id}"
    exit 1
fi

