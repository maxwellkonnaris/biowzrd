#!/bin/bash
#SBATCH --job-name=parallel_rsync
#SBATCH --output=rsync_%j.log
#SBATCH --error=rsync_%j.err
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=ALL

###############################################################################
# This script copies an entire directory (and all its subdirectories) from
# SRC_DIR to DEST_DIR, in parallel, using rsync. It does not delete or modify
# files in the source. At the end, it checks that the destination has the same
# number of files, and logs any missing ones if there's a mismatch.
###############################################################################

set -eu  # Exit on error or undefined variable

# --------------------- Defaults / Usage --------------------
DEFAULT_SRC_DIR="/storage/home/mak6930/scratch/mlscale/fastq_data"
DEFAULT_DEST_DIR="/storage/home/mak6930/silvermanlab/mlscale/fastq_data"

usage() {
cat <<EOF
Usage: $(basename "$0") [options]

Copies an entire directory (with subdirectories) from SRC_DIR to DEST_DIR
using rsync in parallel. Does NOT delete source files.

Options:
  -s <path>  Source directory (default: $DEFAULT_SRC_DIR)
  -d <path>  Destination directory (default: $DEFAULT_DEST_DIR)
  -h         Show this help message and exit

Slurm job parameters (time, CPUs, memory) can be customized in the SBATCH header.
EOF
exit 0
}

SRC_DIR="$DEFAULT_SRC_DIR"
DEST_DIR="$DEFAULT_DEST_DIR"

# ---------------------- Parse CLI args ---------------------
while getopts "s:d:h" opt; do
  case "$opt" in
    s) SRC_DIR="$OPTARG" ;;
    d) DEST_DIR="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))

# ---------------------- Validation -------------------------
if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERROR: Source directory does not exist: $SRC_DIR"
  exit 1
fi

mkdir -p "$DEST_DIR" || {
  echo "ERROR: Unable to create destination directory: $DEST_DIR"
  exit 1
}

# Make a local log file
LOG_FILE="rsync_${SLURM_JOB_ID}.log"
touch "$LOG_FILE" || {
  echo "ERROR: Unable to write to current directory for logs."
  exit 1
}

echo "==================== STARTING JOB $(date) ====================" | tee -a "$LOG_FILE"
echo "Source:      $SRC_DIR" | tee -a "$LOG_FILE"
echo "Destination: $DEST_DIR" | tee -a "$LOG_FILE"

# Count total files in source
echo "Counting files in source..." | tee -a "$LOG_FILE"
EXPECTED_FILE_COUNT=$(find "$SRC_DIR" -type f | wc -l)
echo "Expected file count: $EXPECTED_FILE_COUNT" | tee -a "$LOG_FILE"

# ----------------- Create the chunk script -----------------
# Each parallel job will run this to copy exactly one file, preserving structure.
CHUNK_SCRIPT="per_file_rsync.sh"
cat << 'EOF' > "$CHUNK_SCRIPT"
#!/bin/bash
SRC_FILE="$1"
SRC_DIR="$2"
DEST_DIR="$3"
LOG_FILE="$4"

# Compute relative path
REL_PATH=$(realpath --relative-to="$SRC_DIR" "$SRC_FILE")

# Create the parent directory in the destination
mkdir -p "$DEST_DIR/$(dirname "$REL_PATH")"

# Now rsync that single file
# Using -a to preserve attributes, timestamps, etc.; no deletion on source
# --partial and --append-verify are optional if you want partial file resumes
rsync -a "$SRC_FILE" "$DEST_DIR/$REL_PATH"
RC=$?

if [ $RC -ne 0 ]; then
  echo "[ERROR] rsync failed for file: $SRC_FILE (exit code: $RC)" >> "$LOG_FILE"
  exit $RC
fi

echo "[OK] Copied file: $SRC_FILE" >> "$LOG_FILE"
EOF

chmod +x "$CHUNK_SCRIPT"

# -------------------- Generate file list -------------------
# We'll pass each file to per_file_rsync.sh in parallel.
echo "Building list of all source files..." | tee -a "$LOG_FILE"
FILE_LIST="files_to_copy_${SLURM_JOB_ID}.txt"
find "$SRC_DIR" -type f > "$FILE_LIST"

FILE_COUNT=$(wc -l < "$FILE_LIST")
echo "Found $FILE_COUNT files in the source directory (should match expected)." | tee -a "$LOG_FILE"

# -------------- Parallel copy of all files ---------------
echo "Starting parallel rsync..." | tee -a "$LOG_FILE"

# Try to load parallel if available (depends on your HPC setup)
module load parallel 2>/dev/null || true

if command -v parallel &>/dev/null; then
  # Using parallel
  parallel --no-notice -j "$SLURM_CPUS_PER_TASK" \
    "$PWD/$CHUNK_SCRIPT" {} "$SRC_DIR" "$DEST_DIR" "$PWD/$LOG_FILE" \
    :::: "$FILE_LIST"
else
  # Fallback: sequential
  echo "GNU Parallel not found, copying files sequentially..." | tee -a "$LOG_FILE"
  while read -r file; do
    "$PWD/$CHUNK_SCRIPT" "$file" "$SRC_DIR" "$DEST_DIR" "$PWD/$LOG_FILE"
  done < "$FILE_LIST"
fi

echo "Parallel rsync completed." | tee -a "$LOG_FILE"

# -------------- Verify final file count --------------
echo "Verifying file count in the destination..." | tee -a "$LOG_FILE"
ACTUAL_COUNT=$(find "$DEST_DIR" -type f | wc -l)
echo "Destination file count: $ACTUAL_COUNT" | tee -a "$LOG_FILE"

if [ "$ACTUAL_COUNT" -eq "$EXPECTED_FILE_COUNT" ]; then
  echo "File counts match! Transfer succeeded." | tee -a "$LOG_FILE"
else
  echo "WARNING: File counts do NOT match!" | tee -a "$LOG_FILE"
  echo "  Expected: $EXPECTED_FILE_COUNT" | tee -a "$LOG_FILE"
  echo "  Got:      $ACTUAL_COUNT" | tee -a "$LOG_FILE"

  echo "Listing missing files..." | tee -a "$LOG_FILE"
  # We'll sort them and compare
  comm -23 \
    <(find "$SRC_DIR" -type f | sed "s|$SRC_DIR/||" | sort) \
    <(find "$DEST_DIR" -type f | sed "s|$DEST_DIR/||" | sort) \
    | while read -r missing_rel; do
      echo "Missing: $SRC_DIR/$missing_rel" | tee -a "$LOG_FILE"
    done
fi

echo "==================== JOB FINISHED $(date) ====================" | tee -a "$LOG_FILE"
