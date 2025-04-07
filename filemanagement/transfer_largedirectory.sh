#!/bin/bash
#SBATCH --job-name=parallel_transfer
#SBATCH --output=transfer_%j.log
#SBATCH --error=transfer_%j.err
#SBATCH --time=72:00:00          # Increased time for robustness
#SBATCH --account=open
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G               # Increased memory for large transfers
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

# Enable command tracing for debugging
set -x

# Paths
SRC_DIR="/storage/home/mak6930/scratch/mlscale/fastq_data"
DEST_DIR="/storage/home/mak6930/silvermanlab/mlscale/fastq_data"
STATUS_DIR="./transfer_status"
LOG_FILE="$SLURM_SUBMIT_DIR/transfer_${SLURM_JOB_ID}.log"

# Ensure SLURM_SUBMIT_DIR is set
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}

# Function to log messages with timestamps
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Validate paths
if [ ! -d "$SRC_DIR" ]; then
  log_message "ERROR: Source directory $SRC_DIR does not exist"
  exit 1
fi
if ! mkdir -p "$DEST_DIR" "$STATUS_DIR"; then
  log_message "ERROR: Failed to create directories $DEST_DIR or $STATUS_DIR"
  exit 1
fi

# Generate relative file list
cd "$SRC_DIR" || { log_message "ERROR: Cannot cd to $SRC_DIR"; exit 1; }
find . -type f -not -name "*.tmp" -not -name "*.swp" > "$SLURM_SUBMIT_DIR/filelist.txt"
if [ ! -s "$SLURM_SUBMIT_DIR/filelist.txt" ]; then
  log_message "ERROR: No files found in $SRC_DIR"
  exit 1
fi
cd "$SLURM_SUBMIT_DIR" || { log_message "ERROR: Cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }

# Split file list into chunks
TOTAL_FILES=$(wc -l < filelist.txt)
CHUNK_SIZE=$(( (TOTAL_FILES + SLURM_CPUS_PER_TASK - 1) / SLURM_CPUS_PER_TASK ))
split -l "$CHUNK_SIZE" filelist.txt filelist_part_ || {
  log_message "ERROR: Failed to split filelist.txt"
  exit 1
}

# Create chunk transfer script
cat << 'EOF' > transfer_chunk.sh
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
  if rsync -avh --progress --stats \
      --partial --append-verify --checksum \
      --files-from="$chunk" \
      "$SRC_DIR/" "$DEST_DIR/"; then
    touch "$status_file"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chunk $chunk completed successfully" | tee -a "$LOG_FILE"
    exit 0
  else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Attempt $attempt failed for chunk $chunk" | tee -a "$LOG_FILE"
    [ $attempt -eq 3 ] && {
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Chunk $chunk failed after 3 attempts" | tee -a "$LOG_FILE"
      exit 1
    }
    sleep 10
  fi
done
EOF

chmod +x transfer_chunk.sh

# Run in parallel with no tty prompts
log_message "===== STARTING PARALLEL TRANSFER ====="
parallel --no-notice --plain --eta --retries 3 -j "$SLURM_CPUS_PER_TASK" \
  ./transfer_chunk.sh {} "$SRC_DIR" "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" ::: filelist_part_* || {
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
rsync -avh --dry-run --checksum "$SRC_DIR/" "$DEST_DIR/" | tee "$SLURM_SUBMIT_DIR/validation_${SLURM_JOB_ID}.txt"
if grep -q "would be" "$SLURM_SUBMIT_DIR/validation_${SLURM_JOB_ID}.txt"; then
  log_message "WARNING: Validation found discrepancies"
else
  log_message "Validation passed: All files match"
fi

# Summary using wc instead of ls to avoid alias issues
log_message "===== JOB COMPLETED ====="
COMPLETED_CHUNKS=$(find "$STATUS_DIR" -name "*.done" | wc -l)
TOTAL_CHUNKS=$(find . -name "filelist_part_*" | wc -l)
log_message "Chunks completed: $COMPLETED_CHUNKS / $TOTAL_CHUNKS"

# Optional cleanup (uncomment to enable)
# rm -f filelist.txt filelist_part_* transfer_chunk.sh "$SLURM_SUBMIT_DIR/validation_${SLURM_JOB_ID}.txt"
