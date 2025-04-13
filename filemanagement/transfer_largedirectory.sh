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

################################################################################
# This script transfers files from a source directory to a destination directory
# in parallel using GNU Parallel and rsync with retries. By default, it copies
# all subdirectories, but you can skip subdirectories with `-n`.
#
# HPC usage notes:
#   - Requests 16 CPUs, 64 GB memory, 48 hours by default.
#   - Adjust concurrency in the parallel command if needed (e.g. -j 4).
#   - A final "full directory sync" is done if subdirectories are included,
#     ensuring everything matches even if a chunk missed something.
################################################################################

set -x  # Enable command tracing for debugging

# ---------------------- Default Values ----------------------
DEFAULT_SRC_DIR="/storage/home/mak6930/scratch/mlscale/fastq_data"
DEFAULT_DEST_DIR="/storage/home/mak6930/silvermanlab/mlscale/fastq_data"
INCLUDE_SUBDIRS=true  # Default is "transfer all subdirectories"
PARTIAL_DIR="${TMPDIR:-/tmp}/rsync_partials"  # Where partial files go
mkdir -p "$PARTIAL_DIR"

# ---------------------- Usage Function ----------------------
usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Transfers files/directories from SRC to DEST using rsync and parallel with retries.

Options:
  -s <path>   Source directory (default: $DEFAULT_SRC_DIR)
  -d <path>   Destination directory (default: $DEFAULT_DEST_DIR)
  -n          Do NOT include subdirectories (only transfer top-level files)
  -h          Show this help message and exit

Environment / HPC Info:
  - Slurm job requesting 16 CPUs, 64 GB RAM, 48h.
  - Adjust concurrency by editing parallel -j in the script or changing --cpus-per-task.
  - If subdirectories are included (default), a final deep rsync is done to catch any
    missed files. If -n is used, that deep sync is skipped.

EOF
  exit 0
}

# ------------------- Parse Command-Line Flags -------------------
while getopts "s:d:nh" opt; do
  case "$opt" in
    s) SRC_DIR="$OPTARG" ;;
    d) DEST_DIR="$OPTARG" ;;
    n) INCLUDE_SUBDIRS=false ;;
    h) usage ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))

# If not set by flags, use defaults
: "${SRC_DIR:=$DEFAULT_SRC_DIR}"
: "${DEST_DIR:=$DEFAULT_DEST_DIR}"

# ------------------- Environment Setup -------------------
SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
STATUS_DIR="./transfer_status"
LOG_FILE="$SLURM_SUBMIT_DIR/transfer_${SLURM_JOB_ID}.log"

mkdir -p "$(dirname "$LOG_FILE")" || {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Cannot create log file directory"
  exit 1
}

# ------------------- Logging & Cleanup -------------------
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

cleanup() {
  rm -f "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt" "$SLURM_SUBMIT_DIR/transfer_chunk.sh" 2>/dev/null
  log_message "Cleaned up temporary files in $STATUS_DIR"
}

# Trap common signals (SIGTERM, SIGINT) for HPC job cancellations
trap "log_message 'Caught SIGTERM, cleaning up...' && cleanup && exit 1" SIGTERM
trap "log_message 'Caught SIGINT, cleaning up...' && cleanup && exit 1" SIGINT
trap cleanup EXIT

# ------------------- Validate Paths -------------------
if [ ! -d "$SRC_DIR" ]; then
  log_message "ERROR: Source directory $SRC_DIR does not exist"
  exit 1
fi

if ! mkdir -p "$DEST_DIR" "$STATUS_DIR"; then
  log_message "ERROR: Failed to create directories $DEST_DIR or $STATUS_DIR"
  exit 1
fi

touch "$STATUS_DIR/test_write" && rm "$STATUS_DIR/test_write" || {
  log_message "ERROR: $STATUS_DIR is not writable"
  exit 1
}

# ------------------- Initial Checks -------------------
log_message "===== INITIAL TRANSFER CHECK ====="
src_files=$(find "$SRC_DIR" -type f | wc -l)
src_dirs=$(find "$SRC_DIR" -type d | wc -l)
log_message "Expected to transfer: $src_files files, $src_dirs directories"
log_message "Top-level items:"
while IFS= read -r -d '' item; do
  log_message "  $item"
done < <(find "$SRC_DIR" -maxdepth 1 -not -path "$SRC_DIR" -print0)

# Force reprocessing by removing old status flags
rm -f "$STATUS_DIR"/*.done 2>/dev/null

if [ "$src_files" -eq 0 ] && [ "$src_dirs" -eq 1 ]; then
  log_message "WARNING: Source directory $SRC_DIR appears empty; proceeding anyway."
fi

# Attempt to load parallel if available
module load parallel 2>/dev/null || true

# ------------------- Create Chunk Script -------------------
cat << 'CHUNKEOF' > "$SLURM_SUBMIT_DIR/transfer_chunk.sh"
#!/bin/bash
src_path="$1"
dest_dir="$2"
status_dir="$3"
log_file="$4"
source_dir="$5"
partial_dir="$6"

# Figure out a relative path for status tracking
rel_path=$(realpath --relative-to="$source_dir" "$src_path" 2>/dev/null || basename "$src_path")
status_file="${status_dir}/$(echo "$rel_path" | tr '/' '_').done"

# If it's already transferred, skip
if [ -f "$status_file" ]; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping $src_path (already completed)" | tee -a "$log_file"
  exit 0
fi

# Create parent directories in the destination (for nested paths)
if [ -d "$src_path" ]; then
  mkdir -p "$dest_dir/$(realpath --relative-to="$source_dir" "$src_path")"
else
  parent_subdir=$(dirname "$(realpath --relative-to="$source_dir" "$src_path")")
  mkdir -p "$dest_dir/$parent_subdir"
fi

for attempt in {1..3}; do
  rsync -a --progress --stats \
        --partial --partial-dir="$partial_dir" \
        --append-verify --checksum \
        "$src_path" "$dest_dir/"
  rc=$?

  if [ $rc -ne 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Attempt $attempt failed for $src_path (rsync exit code: $rc)" | tee -a "$log_file"
    if [ $attempt -eq 3 ]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Transfer of $src_path failed after 3 attempts" | tee -a "$log_file"
      exit 1
    fi
    sleep 10
    continue
  fi

  # Success
  touch "$status_file"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Transfer of $src_path completed successfully" | tee -a "$log_file"
  exit 0
done
CHUNKEOF

chmod +x "$SLURM_SUBMIT_DIR/transfer_chunk.sh"

# ------------------- Collect Items to Transfer -------------------
log_message "===== COLLECTING ITEMS TO TRANSFER ====="
chunk_list=()

if $INCLUDE_SUBDIRS; then
  # Include both files and directories at top-level
  while IFS= read -r -d '' item; do
    chunk_list+=("$item")
  done < <(find "$SRC_DIR" -maxdepth 1 -not -path "$SRC_DIR" -print0)
else
  # Only grab top-level files, skip directories
  while IFS= read -r -d '' item; do
    chunk_list+=("$item")
  done < <(find "$SRC_DIR" -maxdepth 1 -type f -print0)
fi

if [ ${#chunk_list[@]} -eq 0 ]; then
  log_message "WARNING: No top-level items to process."
fi

# ------------------- Transfer in Parallel (or Sequential) -------------------
log_message "===== STARTING TRANSFER ====="
if command -v parallel >/dev/null 2>&1; then
  log_message "Using GNU Parallel for multi-core transfers."
  parallel --no-notice --plain --eta --retries 3 -j "$SLURM_CPUS_PER_TASK" \
    "$SLURM_SUBMIT_DIR/transfer_chunk.sh" {} "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" "$SRC_DIR" "$PARTIAL_DIR" ::: "${chunk_list[@]}" || {
    log_message "ERROR: Parallel transfer failed"
    exit 1
  }
else
  log_message "GNU Parallel not found, falling back to sequential transfers."
  for item in "${chunk_list[@]}"; do
    "$SLURM_SUBMIT_DIR/transfer_chunk.sh" "$item" "$DEST_DIR" "$STATUS_DIR" "$LOG_FILE" "$SRC_DIR" "$PARTIAL_DIR" || {
      log_message "ERROR: Transfer of $item failed"
      exit 1
    }
  done
fi

# ------------------- Optional Final Sync & Validation -------------------
if $INCLUDE_SUBDIRS; then
  # If subdirectories are allowed, do a final deep sync to catch any missed items
  log_message "===== FINAL SYNC ====="
  rsync -a --progress --checksum \
        --partial --partial-dir="$PARTIAL_DIR" --append-verify \
        "$SRC_DIR/" "$DEST_DIR/" || {
    log_message "ERROR: Final sync failed"
    exit 1
  }
else
  log_message "Skipping final deep sync (no subdirectories requested)."
fi

# -------------- Final Checks --------------
log_message "===== FINAL TRANSFER CHECK ====="
dest_files=$(find "$DEST_DIR" -type f | wc -l)
dest_dirs=$(find "$DEST_DIR" -type d | wc -l)
log_message "Transferred: $dest_files files, $dest_dirs directories"

if [ "$src_files" -ne "$dest_files" ] || [ "$src_dirs" -ne "$dest_dirs" ]; then
  log_message "ERROR: Mismatch (Expected $src_files files/$src_dirs dirs, got $dest_files files/$dest_dirs dirs)."
  log_message "Checking for missing files..."
  comm -23 <(find "$SRC_DIR" -type f | sort) <(find "$DEST_DIR" -type f | sort) | while read -r missing; do
    log_message "Missing file: $missing"
  done
  log_message "Checking for missing directories..."
  comm -23 <(find "$SRC_DIR" -type d | sort) <(find "$DEST_DIR" -type d | sort) | while read -r missing; do
    log_message "Missing directory: $missing"
  done
  exit 1
fi

# Dry-run validation with checksum
rsync -a --dry-run --checksum "$SRC_DIR/" "$DEST_DIR/" > "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"
if grep -qE "would be|differ" "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt"; then
  log_message "ERROR: Validation found discrepancies"
  cat "$STATUS_DIR/validation_${SLURM_JOB_ID}.txt" >> "$LOG_FILE"
  exit 1
else
  log_message "Validation passed: All files and directories match"
fi

# -------------- Summary --------------
log_message "===== JOB COMPLETED ====="
COMPLETED_ITEMS=$(find "$STATUS_DIR" -name "*.done" | wc -l)
TOTAL_ITEMS=${#chunk_list[@]}
log_message "Top-level items processed: $COMPLETED_ITEMS / $TOTAL_ITEMS"
