#!/bin/bash
#SBATCH --job-name=parallel_transfer
#SBATCH --output=transfer_%j.log
#SBATCH --error=transfer_%j.err
#SBATCH --time=48:00:00
#SBATCH --account=one
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G  # Reduced memory to a more reasonable amount
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

# Enable command tracing for debugging
set -x

# Default paths (used if no arguments provided)
DEFAULT_SRC_DIR="/storage/home/mak6930/scratch/mlscale/fastq_data"
DEFAULT_DEST_DIR="/storage/home/mak6930/silvermanlab/mlscale/fastq_data"

# Set paths from arguments or defaults
SRC_DIR="${1:-$DEFAULT_SRC_DIR}"
DEST_DIR="${2:-$DEFAULT_DEST_DIR}"
STATUS_DIR="./transfer_status"
LOG_FILE="$SLURM_SUBMIT_DIR/transfer_${SLURM_JOB_ID}.log"

# Ensure SLURM_SUBMIT_DIR is set
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}

# Ensure log file directory exists
mkdir -p "$(dirname "$LOG_FILE")" || { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Cannot create log file directory"; exit 1; }

# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Cleanup function for filelists and temp files
cleanup() {
  rm -f "$STATUS_DIR/filelist_part_"* "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt" \
        "$SLURM_SUBMIT_DIR/transfer_chunk.sh" "$STATUS_DIR"/*.done 2>/dev/null
  log_message "Cleaned up temporary files in $STATUS_DIR"
}
trap cleanup EXIT  # Run cleanup on script exit (success, failure, or interrupt)

# Validate paths
if [ ! -d "$SRC_DIR" ]; then
  log_message "ERROR: Source directory $SRC_DIR does not exist"
  exit 1
fi
if ! mkdir -p "$DEST_DIR" "$STATUS_DIR"; then
  log_message "ERROR: Failed to create directories $DEST_DIR or $STATUS_DIR"
  exit 1
fi

# Check if STATUS_DIR is writable
touch "$STATUS_DIR/test_write" && rm "$STATUS_DIR/test_write" || {
  log_message "ERROR: $STATUS_DIR is not writable"
  exit 1
}

# Generate relative file list, stripping ./ prefix
cd "$SRC_DIR" || { log_message "ERROR: Cannot cd to $SRC_DIR"; exit 1; }
find . -type f -not -name "*.tmp" -not -name "*.swp" | sed 's|^\./||' > "$STATUS_DIR/filelist.txt"
if [ ! -s "$STATUS_DIR/filelist.txt" ]; then
  log_message "ERROR: No files found in $SRC_DIR"
  exit 1
fi
cd "$SLURM_SUBMIT_DIR" || { log_message "ERROR: Cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

# Split file list into chunks in STATUS_DIR
TOTAL_FILES=$(wc -l < "$STATUS_DIR/filelist.txt")
CHUNK_SIZE=$(( (TOTAL_FILES + SLURM_CPUS_PER_TASK - 1) / SLURM_CPUS_PER_TASK ))
split -l "$CHUNK_SIZE" "$STATUS_DIR/filelist.txt" "$STATUS_DIR/filelist_part_" || {
  log_message "ERROR: Failed to split filelist.txt"
  exit 1
}
rm -f "$STATUS_DIR/filelist.txt"  # Remove original filelist after splitting

# Create chunk transfer script
cat << 'EOF' > "$SLURM_SUBMIT_DIR/transfer_chunk.sh"
#!/bin/bash
chunk="$1"
SRC_DIR="$2"
DEST_DIR="$3"
STATUS_DIR="$4"
LOG_FILE="$5"

status_file="${STATUS_DIR}/$(basename "$chunk").done"

# Skip if already completed
if [ -f "$status_file" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping $chunk (already completed)" | tee -a "$LOG_FILE"
  exit 0
fi

# Retry logic for robustness
for attempt in {1..3}; do
  rsync -avh --progress --stats \
      --partial --append-verify --checksum \
      --files-from="$chunk" \
      "$SRC_DIR/" "$DEST_DIR/"
  rsync_exit=$?
  if [ $rsync_exit -eq 0 ]; then
    touch "$status_file"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chunk $chunk completed successfully" | tee -a "$LOG_FILE"
    exit 0
  else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Attempt $attempt failed for chunk $chunk (rsync exit code: $rsync_exit)" | tee -a "$LOG_FILE"
    [ $attempt -eq 3 ] && {
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Chunk $chunk failed after 3 attempts" | tee -a "$LOG_FILE"
      exit 1
    }
    sleep 10
  fi
done
EOF

chmod +x "$SLURM_SUBMIT_DIR/transfer_chunk.sh"

# Run in parallel with no tty prompts, skipping empty chunks
log_message "===== STARTING PARALLEL TRANSFER ====="
chunk_list=()
for chunk in "$STATUS_DIR/filelist_part_"*; do
  [ -s "$chunk" ] || { log_message "WARNING: Skipping empty chunk $chunk"; continue; }
  chunk_list+=("$chunk")
done
if [ ${#chunk_list[@]} -eq 0 ]; then
  log_message "ERROR: No valid chunks to process"
  exit 1
fi
parallel --no-notice --plain --eta --retries 3 -j "$SLURM_CPUS_PER_TASK" \
  "$SLURM_SUBMIT_DIR/transfer_chunk.sh" {} "$SRC_DIR" "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" ::: "${chunk_list[@]}" || {
  log_message "ERROR: Parallel transfer failed"
  exit 1
}

# SAFER final sync with verification
log_message "===== FINAL SYNC ====="
rsync -avh --progress --checksum "$SRC_DIR/" "$DEST_DIR/" || {
  log_message "ERROR: Final sync failed"
  exit 1
}

# Validation step
log_message "===== VALIDATING TRANSFER ====="
src_count=$(find "$SRC_DIR" -type f | wc -l)
dest_count=$(find "$DEST_DIR" -type f | wc -l)
if [ "$src_count" -ne "$dest_count" ]; then
  log_message "ERROR: File count mismatch (Source: $src_count, Destination: $dest_count)"
  exit 1
fi
rsync -avh --dry-run --checksum "$SRC_DIR/" "$DEST_DIR/" > "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"
if grep -qE "would be|differ" "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"; then
  log_message "ERROR: Validation found discrepancies"
  exit 1
else
  log_message "Validation passed: All files match"
fi

# Summary using find instead of ls
log_message "===== JOB COMPLETED ====="
COMPLETED_CHUNKS=$(find "$STATUS_DIR" -name "*.done" | wc -l)
TOTAL_CHUNKS=$(find "$STATUS_DIR" -name "filelist_part_*" | wc -l)
log_message "Chunks completed: $COMPLETED_CHUNKS / $TOTAL_CHUNKS"
