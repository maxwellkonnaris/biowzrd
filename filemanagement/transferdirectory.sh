#!/bin/bash
#SBATCH --job-name=parallel_transfer
#SBATCH --output=transfer_%j.log
#SBATCH --error=transfer_%j.err
#SBATCH --time=48:00:00
#SBATCH --account=one
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16        
#SBATCH --mem=16G                   
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

# Source and destination paths
SRC_DIR="/path/to/source_directory"
DEST_DIR="."

# Validate paths
if [ ! -d "$SRC_DIR" ]; then
  echo "ERROR: Source directory $SRC_DIR does not exist" >&2
  exit 1
fi
mkdir -p "$DEST_DIR" || { echo "Failed to create $DEST_DIR"; exit 1; }

# Generate a list of files to transfer (excluding temp files)
find "$SRC_DIR" -type f -not -name "*.tmp" -not -name "*.swp" > filelist.txt

# Split the file list into chunks (1 chunk per CPU core)
split -n l/$SLURM_CPUS_PER_TASK filelist.txt filelist_part_

# Function to transfer a chunk of files
transfer_chunk() {
  chunk="$1"
  rsync -avh --progress --stats \
    --partial --append-verify \
    --files-from="$chunk" \
    "$SRC_DIR/" "$DEST_DIR/"
  echo "Chunk $chunk completed: $(date)"
}
export -f transfer_chunk

# Run transfers in parallel using GNU Parallel
echo "===== STARTING PARALLEL TRANSFER: $(date) ====="
parallel -j $SLURM_CPUS_PER_TASK transfer_chunk ::: filelist_part_*

# Final full sync to catch any missed files
echo "===== FINAL SYNC: $(date) ====="
rsync -avh --delete --progress "$SRC_DIR/" "$DEST_DIR/"

echo "===== JOB COMPLETED: $(date) ====="
