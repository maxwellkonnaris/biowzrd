#!/bin/bash
#SBATCH --job-name=parallel_rsync
#SBATCH --output=rsync_%j.log
#SBATCH --error=rsync_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

################################################################################
# This script copies all files from SRC_DIR to DEST_DIR in parallel with rsync.
# It never deletes the source. Each file's completion is "checkpointed" so if
# the job is stopped and rerun, it will skip already-copied files. At the end,
# it verifies that all files are present. Missing files, if any, get logged.
#
# Key features:
#   - Parallel using GNU Parallel
#   - Checkpointing with per-file .done markers
#   - --partial and --append-verify allow large-file resume
#   - Final mismatch check lists missing files
################################################################################

set -euo pipefail  # strict mode: exit on error, no unset vars, pipeline errors fail

# --------------------- Defaults & Usage --------------------
DEFAULT_SRC_DIR="/storage/home/mak6930/scratch/mlscale/fastq_data"
DEFAULT_DEST_DIR="/storage/home/mak6930/silvermanlab/mlscale/fastq_data"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Copies an entire directory (subdirectories + files) from SRC_DIR to DEST_DIR
in parallel, with robust checkpointing.

Options:
  -s <path>  Source directory (default: $DEFAULT_SRC_DIR)
  -d <path>  Destination directory (default: $DEFAULT_DEST_DIR)
  -h         Show this help message and exit

Slurm job parameters (time, CPUs, memory) can be customized in the SBATCH header.
To resume a partially completed run, just submit the same command again.
EOF
  exit 0
}

SRC_DIR="$DEFAULT_SRC_DIR"
DEST_DIR="$DEFAULT_DEST_DIR"

while getopts "s:d:h" opt; do
  case "$opt" in
    s) SRC_DIR="$OPTARG" ;;
    d) DEST_DIR="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))

# ---------------------- Validate Paths ----------------------
if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERROR: Source directory does not exist: $SRC_DIR"
  exit 1
fi

mkdir -p "$DEST_DIR" || {
  echo "ERROR: Cannot create destination directory: $DEST_DIR"
  exit 1
}

LOG_FILE="rsync_${SLURM_JOB_ID}.log"
touch "$LOG_FILE" || {
  echo "ERROR: Cannot write log file in current directory."
  exit 1
}

STATUS_DIR="./transfer_status_${SLURM_JOB_ID}"
mkdir -p "$STATUS_DIR"

echo "===== JOB START: $(date) =====" | tee -a "$LOG_FILE"
echo "Source:      $SRC_DIR" | tee -a "$LOG_FILE"
echo "Destination: $DEST_DIR" | tee -a "$LOG_FILE"
echo "Status dir:  $STATUS_DIR" | tee -a "$LOG_FILE"

# ----- Count total files in source -----
echo "Counting total files in source..." | tee -a "$LOG_FILE"
EXPECTED_FILE_COUNT=$(find "$SRC_DIR" -type f | wc -l)
echo "Expected file count: $EXPECTED_FILE_COUNT" | tee -a "$LOG_FILE"

# ----------------- Create the chunk script -----------------
# Each parallel job copies exactly 1 file if it's not done yet.
CHUNK_SCRIPT="copy_one_file.sh"
cat << 'EOF' > "$CHUNK_SCRIPT"
#!/bin/bash
SRC_FILE="$1"
SRC_DIR="$2"
DEST_DIR="$3"
STATUS_DIR="$4"
LOG_FILE="$5"

# Derive a relative path
REL_PATH=$(realpath --relative-to="$SRC_DIR" "$SRC_FILE")
DONE_FILE="${STATUS_DIR}/$(echo "$REL_PATH" | tr '/' '_').done"

# If it's already "done", skip
if [[ -f "$DONE_FILE" ]]; then
  echo "[SKIP] Already copied: $SRC_FILE" >> "$LOG_FILE"
  exit 0
fi

# Ensure parent directories exist at the destination
mkdir -p "$DEST_DIR/$(dirname "$REL_PATH")"

# Attempt rsync up to 3 times
for attempt in 1 2 3; do
  rsync -a \
        --partial --append-verify \
        "$SRC_FILE" "$DEST_DIR/$REL_PATH"
  RC=$?

  if [[ $RC -eq 0 ]]; then
    # Success: checkpoint
    touch "$DONE_FILE"
    echo "[OK] Copied file: $SRC_FILE" >> "$LOG_FILE"
    exit 0
  else
    echo "[ERROR] Attempt $attempt failed for: $SRC_FILE (exit code: $RC)" >> "$LOG_FILE"
    if [[ $attempt -eq 3 ]]; then
      echo "[FATAL] Giving up on: $SRC_FILE" >> "$LOG_FILE"
      exit $RC
    fi
    sleep 5
  fi
done
EOF

chmod +x "$CHUNK_SCRIPT"

# ------------- Build a list of ALL files to copy -------------
FILE_LIST="files_to_copy_${SLURM_JOB_ID}.txt"
echo "Generating file list..." | tee -a "$LOG_FILE"
find "$SRC_DIR" -type f > "$FILE_LIST"
TOTAL_FILES=$(wc -l < "$FILE_LIST")
echo "Found $TOTAL_FILES files (should match expected: $EXPECTED_FILE_COUNT)" | tee -a "$LOG_FILE"

# -------------- Parallel copy of all files --------------
echo "Starting parallel rsync..." | tee -a "$LOG_FILE"

module load parallel 2>/dev/null || true
if command -v parallel &>/dev/null; then
  parallel --no-notice -j "$SLURM_CPUS_PER_TASK" \
    "$PWD/$CHUNK_SCRIPT" {} "$SRC_DIR" "$DEST_DIR" "$STATUS_DIR" "$PWD/$LOG_FILE" \
    :::: "$FILE_LIST"
else
  echo "GNU Parallel not found; doing sequential copy." | tee -a "$LOG_FILE"
  while read -r f; do
    "$PWD/$CHUNK_SCRIPT" "$f" "$SRC_DIR" "$DEST_DIR" "$STATUS_DIR" "$PWD/$LOG_FILE"
  done < "$FILE_LIST"
fi

echo "Parallel transfer complete." | tee -a "$LOG_FILE"

# -------------- Verify final file count --------------
echo "Checking file count in destination..." | tee -a "$LOG_FILE"
ACTUAL_FILE_COUNT=$(find "$DEST_DIR" -type f | wc -l)
echo "Destination has $ACTUAL_FILE_COUNT files." | tee -a "$LOG_FILE"

if [[ "$ACTUAL_FILE_COUNT" -eq "$EXPECTED_FILE_COUNT" ]]; then
  echo "All files transferred successfully!" | tee -a "$LOG_FILE"
else
  echo "WARNING: Missing some files (expected $EXPECTED_FILE_COUNT, got $ACTUAL_FILE_COUNT)" | tee -a "$LOG_FILE"
  echo "Listing missing files..." | tee -a "$LOG_FILE"
  comm -23 \
    <(find "$SRC_DIR" -type f | sort) \
    <(find "$DEST_DIR" -type f | sort) \
    | while read -r missing_file; do
      echo "Missing: $missing_file" | tee -a "$LOG_FILE"
    done
fi

echo "===== JOB FINISH: $(date) =====" | tee -a "$LOG_FILE"
