#!/bin/bash
#SBATCH --job-name=process_parallel
#SBATCH --output=process_parallel_%j.out
#SBATCH --error=process_parallel_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

# Default directory for FASTQ files
DEFAULT_DIR="fastq_data/fastq_biologicaldata/"
DEBUG_FILE="debug.log"
DEBUG_LOCK="debug.lock"
COMPLETED_FILE="completed_steps.log"
COMPLETED_LOCK="completed.lock"
FAILED_FILE="failed.log"
FAILED_LOCK="failed.lock"
DEFAULT_WORKERS=8  # Default number of workers

# Parse command-line arguments
while getopts ":i:d:w:" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    d) FASTQ_DIR="$OPTARG" ;;
    w) NUM_WORKERS="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Set defaults
FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_DIR}"
NUM_WORKERS="${NUM_WORKERS:-$DEFAULT_WORKERS}"

# Validate input file and directories
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv> [-d <fastq_directory>] [-w <num_workers>]"
  echo "  <input.tsv>: Tab-separated file with columns: Bioproject, RunAccession, SequencingType"
  echo "  <fastq_directory>: Directory containing FASTQ files (default: $DEFAULT_DIR)"
  echo "  <num_workers>: Number of parallel workers (default: $DEFAULT_WORKERS)"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Input file $INPUT_FILE does not exist."
  exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

# Calculate threads per worker
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
if [[ $THREADS_PER_WORKER -lt 1 ]]; then
  THREADS_PER_WORKER=1
  NUM_WORKERS=$SLURM_CPUS_PER_TASK
fi
echo "Running with $NUM_WORKERS workers, $THREADS_PER_WORKER threads each (total $SLURM_CPUS_PER_TASK CPUs)"

# Logging function with file locking
log_debug() {
  local msg="$1"
  echo "$(date) $msg" >&2
  flock -x "$DEBUG_LOCK" -c "echo \"$(date) $msg\" >> \"$DEBUG_FILE\"" 2>/dev/null || \
    echo "[WARN] Could not log to $DEBUG_FILE: $msg" >&2
}

# Append to file with locking
append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  flock -x "$lock" -c "echo \"$line\" >> \"$file\"" 2>/dev/null || \
    echo "[ERROR] Could not append '$line' to $file" >&2
}

# Validate FASTQ file (size and gzip integrity)
validate_fastq() {
  local file="$1"
  local min_size=100  # Minimum size in bytes
  if [[ ! -f "$file" ]]; then
    log_debug "[validate_fastq] File $file does not exist"
    return 1
  fi
  local size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    log_debug "[validate_fastq] File $file too small: $size bytes"
    return 1
  fi
  if ! gzip -t "$file" 2>/dev/null; then
    log_debug "[validate_fastq] File $file is not a valid gzip"
    return 1
  fi
  return 0
}

# Run command with retry
run_command() {
  local cmd="$1"
  local err_msg="$2"
  local max_attempts=3
  local attempt=1
  local wait_time=2

  while (( attempt <= max_attempts )); do
    if eval "$cmd"; then
      return 0
    else
      log_debug "$err_msg (attempt $attempt/$max_attempts)"
      if (( attempt == max_attempts )); then
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

# Load completed samples into an associative array
declare -A COMPLETED_SAMPLES
load_completed() {
  if [[ -f "$COMPLETED_FILE" ]]; then
    while IFS=":" read -r accession status; do
      if [[ "$status" == "COMPLETE" ]]; then
        COMPLETED_SAMPLES["$accession"]=1
      fi
    done < "$COMPLETED_FILE"
  fi
  log_debug "Loaded ${#COMPLETED_SAMPLES[@]} completed samples from $COMPLETED_FILE"
}

# Process a single sample with checkpointing
process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SAMPLE_TYPE="$3"

  # Check if already completed
  if [[ -n "${COMPLETED_SAMPLES[$RUN_ACCESSION]}" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already fully processed"
    return
  fi

  # Define file paths
  local INPUT_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_1.fastq.gz"
  local PAIRED_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_2.fastq.gz"
  local OUTPUT_DIR="${BIOPROJECT}"
  local LOG_FILE="${OUTPUT_DIR}/processed_files.log"
  local CHECKPOINT_FILE="${OUTPUT_DIR}/checkpoints_${RUN_ACCESSION}.log"
  mkdir -p "$OUTPUT_DIR"

  # Initialize log files
  touch "$LOG_FILE" "$CHECKPOINT_FILE"

  # Validate input files
  if [[ -f "$INPUT_FASTQ" && -f "$PAIRED_FASTQ" ]]; then
    if ! validate_fastq "$INPUT_FASTQ" || ! validate_fastq "$PAIRED_FASTQ"; then
      log_debug "Invalid FASTQ files for $RUN_ACCESSION: $INPUT_FASTQ or $PAIRED_FASTQ"
      append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
      return 1
    fi
    if grep -q "$INPUT_FASTQ" "$LOG_FILE" && grep -q "$PAIRED_FASTQ" "$LOG_FILE"; then
      log_debug "Skipping $RUN_ACCESSION: files already logged"
    else
      append_with_lock "$INPUT_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
      append_with_lock "$PAIRED_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
    fi
  elif [[ -f "$INPUT_FASTQ" ]]; then
    if ! validate_fastq "$INPUT_FASTQ"; then
      log_debug "Invalid FASTQ file for $RUN_ACCESSION: $INPUT_FASTQ"
      append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
      return 1
    fi
    if grep -q "$INPUT_FASTQ" "$LOG_FILE"; then
      log_debug "Skipping $RUN_ACCESSION: file already logged"
    else
      append_with_lock "$INPUT_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
    fi
  else
    log_debug "No valid FASTQ files for $RUN_ACCESSION in $FASTQ_DIR"
    append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
    return 1
  fi

  # Determine resources (for logging)
  local FILE_SIZE_GB=$(du -BG "$INPUT_FASTQ" | cut -f1 | sed 's/G//')
  local CPUS=$THREADS_PER_WORKER
  local MEM TIME
  if (( FILE_SIZE_GB <= 5 )); then
    MEM="8G"; TIME="12:00:00"
  elif (( FILE_SIZE_GB <= 20 )); then
    MEM="16G"; TIME="24:00:00"
  else
    MEM="32G"; TIME="48:00:00"
  fi
  log_debug "Processing $RUN_ACCESSION (Type: $SAMPLE_TYPE, Threads: $CPUS, Mem: $MEM, Size: $FILE_SIZE_GB GB)"

  # Workflow with checkpointing
  if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    if [[ -f "$PAIRED_FASTQ" ]]; then
      local QC1="${OUTPUT_DIR}/qc_${RUN_ACCESSION}_1.fastq.gz"
      local QC2="${OUTPUT_DIR}/qc_${RUN_ACCESSION}_2.fastq.gz"
      if ! grep -q "^${RUN_ACCESSION}:FASTP$" "$CHECKPOINT_FILE"; then
        run_command "conda run -n dada2 fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER" \
          "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
        append_with_lock "${RUN_ACCESSION}:FASTP" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
      fi
      if validate_fastq "$QC1" && validate_fastq "$QC2" && ! grep -q "^${RUN_ACCESSION}:DADA2$" "$CHECKPOINT_FILE"; then
        run_command "conda run -n dada2 Rscript run_dada2.R \"$QC1\" \"$QC2\"" \
          "[dada2] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
        append_with_lock "${RUN_ACCESSION}:DADA2" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
      fi
    else
      local QC="${OUTPUT_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
      if ! grep -q "^${RUN_ACCESSION}:FASTP$" "$CHECKPOINT_FILE"; then
        run_command "conda run -n dada2 fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER" \
          "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
        append_with_lock "${RUN_ACCESSION}:FASTP" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
      fi
      if validate_fastq "$QC" && ! grep -q "^${RUN_ACCESSION}:DADA2$" "$CHECKPOINT_FILE"; then
        run_command "conda run -n dada2 Rscript run_dada2.R \"$QC\"" \
          "[dada2] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
        append_with_lock "${RUN_ACCESSION}:DADA2" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
      fi
    fi
  elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    local QC="${OUTPUT_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
    if ! grep -q "^${RUN_ACCESSION}:FASTP$" "$CHECKPOINT_FILE"; then
      run_command "conda run -n metaphlan fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER" \
        "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
      append_with_lock "${RUN_ACCESSION}:FASTP" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
    fi
    if validate_fastq "$QC" && ! grep -q "^${RUN_ACCESSION}:METAPHLAN$" "$CHECKPOINT_FILE"; then
      run_command "conda run -n metaphlan metaphlan \"$QC\" --input_type fastq --output \"${OUTPUT_DIR}/metaphlan_${RUN_ACCESSION}.txt\" -t $THREADS_PER_WORKER" \
        "[metaphlan] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
      append_with_lock "${RUN_ACCESSION}:METAPHLAN" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
    fi
    if validate_fastq "$QC" && ! grep -q "^${RUN_ACCESSION}:MOTUS$" "$CHECKPOINT_FILE"; then
      run_command "conda run -n motus motus profile -s \"$QC\" -o \"${OUTPUT_DIR}/motus_${RUN_ACCESSION}.txt\" -t $THREADS_PER_WORKER" \
        "[motus] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; return 1; }
      append_with_lock "${RUN_ACCESSION}:MOTUS" "$CHECKPOINT_FILE" "$COMPLETED_LOCK"
    fi
  else
    log_debug "Invalid sample type: $SAMPLE_TYPE for $RUN_ACCESSION"
    append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
    return 1
  fi

  # Mark as complete
  append_with_lock "${RUN_ACCESSION}:COMPLETE" "$COMPLETED_FILE" "$COMPLETED_LOCK"
  log_debug "Finished processing $RUN_ACCESSION"
}

# Final validation and retry of missing samples
final_validation() {
  log_debug "Starting final validation"

  # Load expected accessions from input.tsv
  declare -A EXPECTED_SAMPLES
  while IFS=$'\t' read -r bioproject accession sample_type; do
    EXPECTED_SAMPLES["$accession"]="$bioproject $sample_type"
  done < <(tail -n +2 "$INPUT_FILE")

  # Find actual processed samples
  declare -A ACTUAL_SAMPLES
  for dir in */; do
    if [[ -d "$dir" && "$dir" != "*/" ]]; then
      for file in "$dir"qc_* "$dir"metaphlan_* "$dir"motus_*; do
        if [[ -f "$file" ]]; then
          local accession=$(basename "$file" | sed -E 's/^(qc|metaphlan|motus)_([^_.]+).*/\2/')
          ACTUAL_SAMPLES["$accession"]="$dir"
        fi
      done
    fi
  done

  # Identify missing samples
  local missing=()
  for accession in "${!EXPECTED_SAMPLES[@]}"; do
    if [[ -z "${ACTUAL_SAMPLES[$accession]}" || -z "${COMPLETED_SAMPLES[$accession]}" ]]; then
      missing+=("$accession")
      log_debug "Missing or incomplete sample: $accession"
    fi
  done

  # Retry missing samples
  if [[ ${#missing[@]} -gt 0 ]]; then
    log_debug "Retrying ${#missing[@]} missing samples"
    printf '%s\n' "${missing[@]}" | parallel --jobs "$NUM_WORKERS" process_sample "${EXPECTED_SAMPLES[{1}]}" {1}
  else
    log_debug "No missing samples to retry"
  fi

  # Log final status
  log_debug "Final validation complete. Expected: ${#EXPECTED_SAMPLES[@]}, Actual: ${#ACTUAL_SAMPLES[@]}, Completed: ${#COMPLETED_SAMPLES[@]}"
}

# Export functions and variables for GNU Parallel
export -f process_sample log_debug append_with_lock validate_fastq run_command final_validation
export FASTQ_DIR DEBUG_FILE DEBUG_LOCK COMPLETED_FILE COMPLETED_LOCK FAILED_FILE FAILED_LOCK THREADS_PER_WORKER

# Initialize lock files
touch "$DEBUG_LOCK" "$COMPLETED_LOCK" "$FAILED_LOCK"

# Load completed samples once
load_completed

# Process input file with GNU Parallel
log_debug "Starting initial processing with $NUM_WORKERS workers"
tail -n +2 "$INPUT_FILE" | parallel --colsep '\t' --jobs "$NUM_WORKERS" process_sample {1} {2} {3}

# Run final validation
final_validation

log_debug "All processing and validation complete."
echo "All processing and validation complete."
