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

# Paths
SRC_DIR="/storage/home/mak6930/scratch/SRA/fastq_data"
DEST_DIR="/storage/home/mak6930/silvermanlab/datarepo/16s"
STATUS_DIR="./transfer_status"

# Validate paths
if [ ! -d "$SRC_DIR" ]; then
  echo "ERROR: Source directory $SRC_DIR does not exist" >&2
  exit 1
fi
mkdir -p "$DEST_DIR" "$STATUS_DIR"

# Generate relative file list
cd "$SRC_DIR" || exit 1
find . -type f -not -name "*.tmp" -not -name "*.swp" > "$SLURM_SUBMIT_DIR/filelist.txt"
cd "$SLURM_SUBMIT_DIR"

# Split file list into chunks
if split --version | grep -q 'GNU coreutils'; then
  split -n l/$SLURM_CPUS_PER_TASK filelist.txt filelist_part_
else
  split -l $(( $(wc -l < filelist.txt) / $SLURM_CPUS_PER_TASK + 1 )) filelist.txt filelist_part_
fi

# Create chunk transfer script
cat << 'EOF' > transfer_chunk.sh
#!/bin/bash
chunk="$1"
SRC_DIR="$2"
DEST_DIR="$3"
STATUS_DIR="$4"

status_file="${STATUS_DIR}/$(basename "$chunk").done"

# Skip if already completed
if [ -f "$status_file" ]; then
  echo "Skipping $chunk (already completed)"
  exit 0
fi

# Perform transfer
if rsync -avh --progress --stats \
    --partial --append-verify \
    --files-from="$chunk" \
    "$SRC_DIR/" "$DEST_DIR/"; then
    touch "$status_file"
    echo "Chunk $chunk completed successfully: $(date)"
else
    echo "ERROR: Chunk $chunk failed: $(date)" >&2
    exit 1
fi
EOF

chmod +x transfer_chunk.sh

# Run in parallel with no tty prompts
echo "===== STARTING PARALLEL TRANSFER: $(date) ====="
parallel --no-notice --plain --eta --retries 3 -j $SLURM_CPUS_PER_TASK \
  ./transfer_chunk.sh {} "$SRC_DIR" "$DEST_DIR" "$STATUS_DIR" ::: filelist_part_*

# Final sync for safety
echo "===== FINAL SYNC: $(date) ====="
rsync -avh --delete --progress "$SRC_DIR/" "$DEST_DIR/"

# Summary
echo "===== JOB COMPLETED: $(date) ====="
echo "Chunks completed: $(ls $STATUS_DIR | wc -l) / $(ls filelist_part_* | wc -l)"

# Optional cleanup (uncomment to enable)
# rm -f filelist.txt filelist_part_* transfer_chunk.sh
