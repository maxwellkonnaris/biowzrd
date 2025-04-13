#!/bin/bash
#SBATCH --job-name=parallel_transfer
#SBATCH --output=transfer_%j.log
#SBATCH --error=transfer_%j.err
#SBATCH --time=48:00:00
#SBATCH --account=one
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

# Enable command tracing for debugging
set -x

# Default paths
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

# Cleanup function
cleanup() {
  rm -f "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt" "$SLURM_SUBMIT_DIR/transfer_chunk.sh" 2>/dev/null
  log_message "Cleaned up temporary files in $STATUS_DIR"
}
trap cleanup EXIT

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

# Initial transfer check
log_message "===== INITIAL TRANSFER CHECK ====="
src_files=$(find "$SRC_DIR" -type f | wc -l)
src_dirs=$(find "$SRC_DIR" -type d | wc -l)
log_message "Expected to transfer: $src_files files, $src_dirs directories"
log_message "Top-level items:"
while IFS= read -r -d '' item; do
  log_message "  $item"
done < <(find "$SRC_DIR" -maxdepth 1 -not -path "$SRC_DIR" -print0)

# Clean up previous .done files to force reprocessing
rm -f "$STATUS_DIR"/*.done 2>/dev/null

# Check if source directory is empty
if [ "$src_files" -eq 0 ] && [ "$src_dirs" -eq 1 ]; then
  log_message "WARNING: Source directory $SRC_DIR is empty, but proceeding to ensure directory structure"
fi

# Try to load parallel module
module load parallel 2>/dev/null || true

# Create chunk transfer script
cat << 'EOF' > "$SLURM_SUBMIT_DIR/transfer_chunk.sh"
#!/bin/bash
src_path="$1"
dest_dir="$2"
status_dir="$3"
log_file="$4"

# Extract relative path for status file
rel_path=$(realpath --relative-to="$SRC_DIR" "$src_path" 2>/dev/null || basename "$src_path")
status_file="${status_dir}/$(echo "$rel_path" | tr '/' '_').done"

# Skip if already completed
if [ -f "$status_file" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping $src_path (already completed)" | tee -a "$log_file"
  exit 0
fi

# Retry logic for robustness
for attempt in {1..3}; do
  rsync -a --progress --stats \
      --partial --append-verify --checksum \
      --timeout=300 \
      "$src_path" "$dest_dir/" || {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Attempt $attempt failed for $src_path (rsync exit code: $?)" | tee -a "$log_file"
    [ $attempt -eq 3 ] && {
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Transfer of $src_path failed after 3 attempts" | tee -a "$log_file"
      exit 1
    }
    sleep 10
    continue
  }
  touch "$status_file"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transfer of $src_path completed successfully" | tee -a "$log_file"
  exit 0
done
EOF

chmod +x "$SLURM_SUBMIT_DIR/transfer_chunk.sh"

# Collect top-level items
log_message "===== COLLECTING ITEMS TO TRANSFER ====="
chunk_list=()
while IFS= read -r -d '' item; do
  chunk_list+=("$item")
done < <(find "$SRC_DIR" -maxdepth 1 -not -path "$SRC_DIR" -print0)
if [ ${#chunk_list[@]} -eq 0 ]; then
  log_message "WARNING: No items to process, but directory structure will be created"
fi

# Transfer items
if command -v parallel >/dev/null 2>&1; then
  log_message "===== STARTING PARALLEL TRANSFER ====="
  parallel --no-notice --plain --eta --retries 3 -j "$SLURM_CPUS_PER_TASK" \
    "$SLURM_SUBMIT_DIR/transfer_chunk.sh" {} "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" ::: "${chunk_list[@]}" || {
    log_message "ERROR: Parallel transfer failed"
    exit 1
  }
else
  log_message "===== STARTING SEQUENTIAL TRANSFER (parallel not found) ====="
  for item in "${chunk_list[@]}"; do
    "$SLURM_SUBMIT_DIR/transfer_chunk.sh" "$item" "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" || {
      log_message "ERROR: Transfer of $item failed"
      exit 1
    }
  done
fi

# Final sync
log_message "===== FINAL SYNC ====="
rsync -a --progress --checksum --timeout=300 "$SRC_DIR/" "$DEST_DIR/" || {
  log_message "ERROR: Final sync failed"
  exit 1
}

# Final validation and comparison
log_message "===== FINAL TRANSFER CHECK ====="
dest_files=$(find "$DEST_DIR" -type f | wc -l)
dest_dirs=$(find "$DEST_DIR" -type d | wc -l)
log_message "Transferred: $dest_files files, $dest_dirs directories"
if [ "$src_files" -ne "$dest_files" ] || [ "$src_dirs" -ne "$dest_dirs" ]; then
  log_message "ERROR: Mismatch detected (Expected: $src_files files, $src_dirs dirs; Got: $dest_files files, $dest_dirs dirs)"
  # Log missing files
  log_message "Checking for missing files..."
  comm -23 <(find "$SRC_DIR" -type f | sort) <(find "$DEST_DIR" -type f | sort) | while read -r missing; do
    log_message "Missing file: $missing"
  done
  # Log missing directories
  log_message "Checking for missing directories..."
  comm -23 <(find "$SRC_DIR" -type d | sort) <(find "$DEST_DIR" -type d | sort) | while read -r missing; do
    log_message "Missing directory: $missing"
  done
  exit 1
fi
rsync -a --dry-run --checksum "$SRC_DIR/" "$DEST_DIR/" > "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"
if grep -qE "would be|differ" "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"; then
  log_message "ERROR: Validation found discrepancies"
  cat "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt" >> "$LOG_FILE"
  exit 1
else
  log_message "Validation passed: All files and directories match"
fi

# Summary
log_message "===== JOB COMPLETED ====="
COMPLETED_ITEMS=$(find "$STATUS_DIR" -name "*.done" | wc -l)
TOTAL_ITEMS=${#chunk_list[@]}
log_message "Top-level items processed: $COMPLETED_ITEMS / $TOTAL_ITEMS"
