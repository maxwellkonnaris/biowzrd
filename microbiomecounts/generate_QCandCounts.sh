#!/bin/bash
export SHELL=/bin/bash

#SBATCH --job-name=counts
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --account=open
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

# Purpose: SLURM script for processing 16S and metagenomic FASTQ files 
#          using fastp (optional), DADA2, MetaPhlAn, and mOTUs.
# Dependencies: micromamba, fastp, R (DADA2), MetaPhlAn, mOTUs, GNU parallel.
# Tip: Use --debug with a small dataset to test the pipeline.

########################################
# Default Variables & Directories
########################################

DEFAULT_DIR="fastq_data/fastq_biologicaldata"
LOCK_DIR="$PWD/locks"
mkdir -p "$LOCK_DIR"

touch "$LOCK_DIR/debug.lock" "$LOCK_DIR/completed.lock" "$LOCK_DIR/failed.lock" "$LOCK_DIR/input.lock"
DEBUG_FILE="debug.log"
DEBUG_LOCK="$LOCK_DIR/debug.lock"
COMPLETED_FILE="completed_steps.log"
COMPLETED_LOCK="$LOCK_DIR/completed.lock"
FAILED_FILE="failed.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
INPUT_LOCK="$LOCK_DIR/input.lock"

DEFAULT_WORKERS=4
DEFAULT_MOTUS_TAX_LEVEL="mOTU"
LOG_LEVEL="INFO"
QC_ENABLED="false"
TMP_BASE="$PWD"
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"

DEFAULT_DADA2_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/dada2"
DEFAULT_MOTUS_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/motus"
DEFAULT_MPA_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa"

INPUT_HEADER=""
declare -a TMP_DIRS

########################################
# Cleanup
########################################
cleanup() {
  for dir in "${TMP_DIRS[@]}"; do
    [[ -d "$dir" ]] && rm -rf "$dir"
  done
  rm -rf "$LOCK_DIR"
}
trap cleanup EXIT

########################################
# Logging Functions
########################################

log_debug() {
  [[ "$LOG_LEVEL" == "DEBUG" ]] || return
  local msg="$1"
  echo "$(date) [DEBUG] $msg" >&2
  {
    flock -x "$DEBUG_LOCK" || true
    echo "$(date) [DEBUG] $msg" >> "$DEBUG_FILE"
  } 2>/dev/null
}

log_info() {
  local msg="$1"
  echo "$(date) [INFO] $msg" | tee -a "$DEBUG_FILE"
}

########################################
# CSV / File Locking Helpers
########################################

append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  {
    flock -x "$lock" || true
    echo "$line" >> "$file"
  } 2>/dev/null
}

update_input_csv() {
  local content="$1"
  local file="$INPUT_FILE"
  local lock="$INPUT_LOCK"
  {
    flock -x "$lock" || true
    printf "%s" "$content" > "$file"
  } 2>/dev/null
}

########################################
# Command-Line Argument Parsing
########################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_FILE="$2"; shift 2 ;;
    -d) FASTQ_DIR="$2"; shift 2 ;;
    -w) NUM_WORKERS="$2"; shift 2 ;;
    -k) MOTUS_TAX_LEVEL="$2"; shift 2 ;;
    -a) DADA2_ENV="$2"; shift 2 ;;
    -m) MOTUS_ENV="$2"; shift 2 ;;
    -p) MPA_ENV="$2"; shift 2 ;;
    --qc) QC_ENABLED="true"; shift ;;
    --debug) LOG_LEVEL="DEBUG"; QC_ENABLED="true"; shift ;;
    -*) echo "Invalid option: $1" >&2; exit 1 ;;
    *) echo "Unexpected argument: $1" >&2; exit 1 ;;
  esac
done

FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_DIR}"
NUM_WORKERS="${NUM_WORKERS:-$DEFAULT_WORKERS}"
MOTUS_TAX_LEVEL="${MOTUS_TAX_LEVEL:-$DEFAULT_MOTUS_TAX_LEVEL}"
DADA2_ENV="${DADA2_ENV:-$DEFAULT_DADA2_ENV}"
MOTUS_ENV="${MOTUS_ENV:-$DEFAULT_MOTUS_ENV}"
MPA_ENV="${MPA_ENV:-$DEFAULT_MPA_ENV}"

if command -v realpath &>/dev/null; then
  FASTQ_DIR="$(realpath "$FASTQ_DIR")"
fi

########################################
# Check micromamba, environments
########################################
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found"
  exit 1
fi

check_env() {
  local env="$1"
  local env_name="$2"
  local env_basename

  if [[ -d "$env" ]]; then
    env_basename=$(basename "$env")
    eval "${env_name}_NAME=$env_basename"
  elif micromamba env list | grep -qE "^[[:space:]]*$env[[:space:]]"; then
    eval "${env_name}_NAME=$env"
  else
    echo "Error: micromamba environment $env_name ($env) not found"
    exit 1
  fi
}

check_env "$DADA2_ENV" "dada2"
check_env "$MOTUS_ENV" "motus"
check_env "$MPA_ENV"   "mpa"

DADA2_ENV_NAME="${dada2_NAME}"
MOTUS_ENV_NAME="${motus_NAME}"
MPA_ENV_NAME="${mpa_NAME}"

########################################
# Validate input file presence
########################################
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv|input.csv|input.txt> [-d <fastq_directory>] [-w <num_workers>] [-k <motus_tax_level>] [--qc] [--debug] ..."
  echo "Tip: Use --debug with a small dataset to test the pipeline."
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file $INPUT_FILE does not exist."
  exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
  echo "Error: Temporary directory base $TMP_BASE does not exist or is not writable."
  exit 1
fi

########################################
# Detect input file format
########################################
validate_input_file() {
  local file="$1"
  local first_line
  first_line=$(head -n 1 "$file")
  
  if [[ "$first_line" == *"Bioproject,RunAccession,SequencingType"* ]]; then
    DELIMITER=','
    EXPECTED_HEADER="Bioproject,RunAccession,SequencingType"
  elif [[ "$first_line" == $'Bioproject\tRunAccession\tSequencingType'* ]]; then
    DELIMITER=$'\t'
    EXPECTED_HEADER="Bioproject\tRunAccession\tSequencingType"
  else
    echo "Error: Input file $file has unrecognized format or header."
    exit 1
  fi

  INPUT_HEADER="$first_line"
  if [[ "$first_line" != "$EXPECTED_HEADER"* ]]; then
    echo "Error: Input file $file has incorrect header. Expected at least: $EXPECTED_HEADER"
    exit 1
  fi
}

validate_input_file "$INPUT_FILE"

########################################
# Check RDP DB if 16S present
########################################
check_rdp_database() {
  local has_16s=false
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    if [[ "$sample_type" == "16S" ]]; then
      has_16s=true
      break
    fi
  done < <(tail -n +2 "$INPUT_FILE")
  
  if [[ "$has_16s" == "true" ]]; then
    if [[ ! -f "$RDP_DATABASE" || ! -r "$RDP_DATABASE" ]]; then
      echo "Error: RDP database $RDP_DATABASE not found or not readable."
      exit 1
    fi
    log_info "RDP database $RDP_DATABASE found for 16S samples."
  fi
}
check_rdp_database

########################################
# SLURM / CPU / Worker Setup
########################################
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
if [[ -z "$NUM_WORKERS" || "$NUM_WORKERS" -le 0 ]]; then
  NUM_WORKERS=$DEFAULT_WORKERS
fi
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
(( THREADS_PER_WORKER < 1 )) && THREADS_PER_WORKER=1

log_info "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK, NUM_WORKERS=$NUM_WORKERS, THREADS_PER_WORKER=$THREADS_PER_WORKER"
log_info "Running with $NUM_WORKERS workers, $THREADS_PER_WORKER threads each"

########################################
# FASTQ Validation
########################################
validate_fastq() {
  local file="$1"
  local min_size=100
  if [[ ! -f "$file" ]]; then
    echo "File does not exist: $file"
    return 1
  fi

  local size
  size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    echo "File too small: $file ($size bytes)"
    return 1
  fi

  if [[ "$file" == *.gz ]]; then
    gzip -t "$file" 2>/dev/null
    if [[ $? -ne 0 ]]; then
      echo "Invalid gzip: $file"
      return 1
    fi
  fi
  return 0
}

########################################
# Command Runners (with retries)
########################################
run_command() {
  local cmd="$1"
  local msg="$2"
  local input_file="$3"
  local output_file="$4"

  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=7200
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command (attempt $attempt): $cmd"
    error_output=$(timeout $timeout bash -c "$cmd" 2>&1)
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

run_command_with_output() {
  local cmd="$1"
  local msg="$2"
  local output_file="$3"
  local input_file="$4"

  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=7200
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command with output (attempt $attempt): $cmd > $output_file"
    error_output=$(timeout $timeout bash -c "$cmd > \"$output_file\" 2>&1")
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

########################################
# Update Input CSV with FASTQ paths
########################################
#!/bin/bash
export SHELL=/bin/bash

#SBATCH --job-name=counts
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --account=open
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --shell=/bin/bash
#SBATCH --open-mode=append

set -o pipefail

# Purpose: SLURM script for processing 16S and metagenomic FASTQ files 
#          using fastp (optional), DADA2, MetaPhlAn, and mOTUs.
# Dependencies: micromamba, fastp, R (DADA2), MetaPhlAn, mOTUs, GNU parallel.
# Tip: Use --debug with a small dataset to test the pipeline.

########################################
# Default Variables & Directories
########################################

DEFAULT_DIR="fastq_data/fastq_biologicaldata"
LOCK_DIR="$PWD/locks"
if ! mkdir -p "$LOCK_DIR" 2>/dev/null; then
  echo "Error: Failed to create lock directory $LOCK_DIR"
  exit 1
fi

touch "$LOCK_DIR/failed.lock" "$LOCK_DIR/input.lock" || {
  echo "Error: Failed to create lock files in $LOCK_DIR"
  exit 1
}
FAILED_FILE="failed.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
INPUT_LOCK="$LOCK_DIR/input.lock"

DEFAULT_WORKERS=4
DEFAULT_MOTUS_TAX_LEVEL="mOTU"
LOG_LEVEL="INFO"
QC_ENABLED="false"
TMP_BASE="$PWD"
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"
METAPHLAN_DB="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa/lib/python3.7/site-packages/metaphlan/metaphlan_databases"

DEFAULT_DADA2_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/dada2"
DEFAULT_MOTUS_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/motus"
DEFAULT_MPA_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa"

INPUT_HEADER=""
declare -a TMP_DIRS

########################################
# Cleanup
########################################
cleanup() {
  for dir in "${TMP_DIRS[@]}"; do
    [[ -d "$dir" ]] && rm -rf "$dir"
  done
}
trap cleanup EXIT

########################################
# Logging Functions
########################################

log_debug() {
  [[ "$LOG_LEVEL" == "DEBUG" ]] || return
  local msg="$1"
  echo "$(date) [DEBUG] $msg" >&2
}

log_info() {
  local msg="$1"
  echo "$(date) [INFO] $msg"
}

########################################
# CSV / File Locking Helpers
########################################

append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  {
    flock -x 200
    echo "$line" >> "$file"
  } 200>"$lock" 2>/dev/null
}

update_input_csv() {
  local content="$1"
  local file="$INPUT_FILE"
  local lock="$INPUT_LOCK"
  {
    flock -x 200
    printf "%s\n" "$content" > "$file"  # Ensure trailing newline
  } 200>"$lock" 2>/dev/null
}

########################################
# Command-Line Argument Parsing
########################################
PROCESS_SAMPLE=0
PROCESS_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_FILE="$2"; shift 2 ;;
    -d) FASTQ_DIR="$2"; shift 2 ;;
    -w) NUM_WORKERS="$2"; shift 2 ;;
    -k) MOTUS_TAX_LEVEL="$2"; shift 2 ;;
    -a) DADA2_ENV="$2"; shift 2 ;;
    -m) MOTUS_ENV="$2"; shift 2 ;;
    -p) MPA_ENV="$2"; shift 2 ;;
    --qc) QC_ENABLED="true"; shift ;;
    --debug) LOG_LEVEL="DEBUG"; shift ;;
    --process-sample)
      PROCESS_SAMPLE=1
      shift
      PROCESS_ARGS=("$@")
      break
      ;;
    -*) echo "Invalid option: $1" >&2; exit 1 ;;
    *) echo "Unexpected argument: $1" >&2; exit 1 ;;
  esac
done

if [[ $PROCESS_SAMPLE -eq 1 ]]; then
  process_sample "${PROCESS_ARGS[@]}"
  exit $?
fi

FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_DIR}"
NUM_WORKERS="${NUM_WORKERS:-$DEFAULT_WORKERS}"
MOTUS_TAX_LEVEL="${MOTUS_TAX_LEVEL:-$DEFAULT_MOTUS_TAX_LEVEL}"
DADA2_ENV="${DADA2_ENV:-$DEFAULT_DADA2_ENV}"
MOTUS_ENV="${MOTUS_ENV:-$DEFAULT_MOTUS_ENV}"
MPA_ENV="${MPA_ENV:-$DEFAULT_MPA_ENV}"
QC_ENABLED="${QC_ENABLED:-false}"

if command -v realpath &>/dev/null; then
  FASTQ_DIR="$(realpath "$FASTQ_DIR")"
fi

########################################
# Check micromamba, environments
########################################
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found"
  exit 1
fi

check_env() {
  local env="$1"
  local env_name="$2"
  local env_basename

  if [[ -d "$env" ]]; then
    env_basename=$(basename "$env")
    eval "${env_name}_NAME=$env_basename"
  elif micromamba env list | grep -qE "^[[:space:]]*$env[[:space:]]"; then
    eval "${env_name}_NAME=$env"
  else
    echo "Error: micromamba environment $env_name ($env) not found"
    exit 1
  fi
}

check_env "$DADA2_ENV" "dada2"
check_env "$MOTUS_ENV" "motus"
check_env "$MPA_ENV"   "mpa"

DADA2_ENV_NAME="${dada2_NAME}"
MOTUS_ENV_NAME="${motus_NAME}"
MPA_ENV_NAME="${mpa_NAME}"

########################################
# Check MetaPhlAn database
########################################
check_metaphlan_db() {
  if [[ -z "$METAPHLAN_DB" ]]; then
    echo "Error: METAPHLAN_DB is not set."
    exit 1
  fi

  log_debug "Checking MetaPhlAn database at $METAPHLAN_DB"
  if ! micromamba run -n "$MPA_ENV_NAME" ls -A "$METAPHLAN_DB"/*.bt2l >/dev/null 2>&1; then
    log_info "MetaPhlAn database not found or empty at $METAPHLAN_DB. Installing..."
    micromamba run -n "$MPA_ENV_NAME" metaphlan --install --bowtie2db "$METAPHLAN_DB" || {
      echo "Error: Failed to install MetaPhlAn database at $METAPHLAN_DB"
      exit 1
    }
    if ! micromamba run -n "$MPA_ENV_NAME" ls -A "$METAPHLAN_DB"/*.bt2l >/dev/null 2>&1; then
      echo "Error: MetaPhlAn database still not found or empty at $METAPHLAN_DB after installation"
      exit 1
    fi
  fi
}
check_metaphlan_db

########################################
# Validate input file presence
########################################
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv|input.csv|input.txt> [-d <fastq_directory>] [-w <num_workers>] [-k <motus_tax_level>] [--qc] [--debug]"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file $INPUT_FILE does not exist."
  exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
  echo "Error: Directory base $TMP_BASE does not exist or is not writable."
  exit 1
fi

########################################
# Detect input file format
########################################
validate_input_file() {
  local file="$1"
  local first_line
  first_line=$(head -n 1 "$file")
  
  if [[ "$first_line" == *"Bioproject,RunAccession,SequencingType"* ]]; then
    DELIMITER=','
    EXPECTED_HEADER="Bioproject,RunAccession,SequencingType"
  elif [[ "$first_line" == $'Bioproject\tRunAccession\tSequencingType'* ]]; then
    DELIMITER=$'\t'
    EXPECTED_HEADER="Bioproject\tRunAccession\tSequencingType"
  else
    echo "Error: Input file $file has unrecognized format or header."
    exit 1
  fi

  INPUT_HEADER="$first_line"
  if [[ "$first_line" != "$EXPECTED_HEADER"* ]]; then
    echo "Error: Input file $file has incorrect header. Expected at least: $EXPECTED_HEADER"
    exit 1
  fi
}

validate_input_file "$INPUT_FILE"

########################################
# Check RDP DB if 16S present
########################################
check_rdp_database() {
  # NEW: Use temporary file instead of process substitution to avoid syntax errors
  local has_16s=false
  if [[ -z "$INPUT_FILE" || ! -f "$INPUT_FILE" ]]; then
    echo "Error: INPUT_FILE is unset or does not exist: $INPUT_FILE"
    exit 1
  fi
  local temp_file
  temp_file=$(mktemp)
  local error_log
  error_log=$(mktemp)
  # NEW: Redirect errors to a log file for debugging
  tail -n +2 "$INPUT_FILE" > "$temp_file" 2> "$error_log"
  if [[ -s "$error_log" ]]; then
    echo "Error: Failed to read $INPUT_FILE: $(cat "$error_log")"
    rm -f "$temp_file" "$error_log"
    exit 1
  fi
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest || [[ -n "$bioproject" ]]; do
    if [[ "$sample_type" == "16S" ]]; then
      has_16s=true
      break
    fi
  done < "$temp_file"
  rm -f "$temp_file" "$error_log"
  if [[ "$has_16s" == "true" ]]; then
    if [[ ! -f "$RDP_DATABASE" || ! -r "$RDP_DATABASE" ]]; then
      echo "Error: RDP database $RDP_DATABASE not found or not readable."
      exit 1
    fi
    log_info "RDP database $RDP_DATABASE found for 16S samples."
  fi
}
# NEW: Explicitly check exit status
check_rdp_database
if [[ $? -ne 0 ]]; then
  echo "Error: check_rdp_database failed"
  exit 1
fi

########################################
# SLURM / CPU / Worker Setup
########################################
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
if [[ -z "$NUM_WORKERS" || "$NUM_WORKERS" -le 0 ]]; then
  NUM_WORKERS=$DEFAULT_WORKERS
fi
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
(( THREADS_PER_WORKER < 1 )) && THREADS_PER_WORKER=1

log_info "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK, NUM_WORKERS=$NUM_WORKERS, THREADS_PER_WORKER=$THREADS_PER_WORKER"

########################################
# FASTQ Validation
########################################
validate_fastq() {
  local file="$1"
  local min_size=100
  if [[ ! -f "$file" ]]; then
    echo "File does not exist: $file"
    return 1
  fi

  local size
  size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    echo "File too small: $file ($size bytes)"
    return 1
  fi

  if [[ "$file" == *.gz ]]; then
    gzip -t "$file" 2>/dev/null
    if [[ $? -ne 0 ]]; then
      echo "Invalid gzip: $file"
      return 1
    fi
  fi
  return 0
}

########################################
# Command Runners (with retries)
########################################
run_command() {
  local cmd="$1"
  local msg="$2"
  local input_file="$3"
  local output_file="$4"

  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=14400
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command (attempt $attempt): $cmd"
    error_output=$(timeout $timeout bash -c "$cmd" 2>&1)
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

run_command_with_output() {
  local cmd="$1"
  local msg="$2"
  local output_file="$3"
  local input_file="$4"

  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=14400
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command with output (attempt $attempt): $cmd > $output_file"
    error_output=$(timeout $timeout bash -c "$cmd > \"$output_file\" 2>&1")
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

########################################
# Update Input CSV with FASTQ paths
########################################
find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq.gz" > fastq_index.txt

# Modified function to use index
update_input_with_fastq_paths() {
  log_info "Checking if input file needs FASTQ path updates"

  local header
  header=$(head -n 1 "$INPUT_FILE")

  # Check if update is needed
  local needs_update=false
  if [[ "$header" == *",Fastq1,Fastq2"* ]]; then
    local all_valid=true
    while IFS="$DELIMITER" read -r b a t f1 f2 _ || [[ -n "$b" ]]; do
      [[ -z "$b" || "$b" == "Bioproject" ]] && continue
      if [[ -z "$f1" && -z "$f2" ]] || [[ -n "$f1" && ! -f "$f1" ]] || [[ -n "$f2" && ! -f "$f2" ]]; then
        all_valid=false
        break
      fi
    done < "$INPUT_FILE"
    [[ "$all_valid" == "true" ]] && { log_info "FASTQ columns valid."; return; }
    needs_update=true
  else
    header="${header},Fastq1,Fastq2"
    needs_update=true
  fi

  if [[ "$needs_update" == "true" ]]; then
    # Check for index file
    if [[ ! -f "fastq_index.txt" ]]; then
      log_info "Error: fastq_index.txt not found. Run: find \"$FASTQ_DIR\" -maxdepth 1 -type f -name \"*.fastq.gz\" > fastq_index.txt"
      exit 1
    fi

    # Update CSV
    local temp_file
    temp_file=$(mktemp)
    echo "$header" > "$temp_file"
    local missing_count=0
    while IFS="$DELIMITER" read -r b a t f1 f2 rest || [[ -n "$b" ]]; do
      [[ -z "$b" || "$b" == "Bioproject" ]] && continue
      local new_fastq1=""
      local new_fastq2=""
      # Search index for accession
      local matches
      mapfile -t matches < <(grep -F "${a}" fastq_index.txt | grep -E "${a}.*\.fastq\.gz$")
      case ${#matches[@]} in
        0)
          ((missing_count++))
          ;;
        1)
          new_fastq1="${matches[0]}"
          ;;
        2)
          new_fastq1="${matches[0]}"  # Assume first is _1 or single-end
          new_fastq2="${matches[1]}"  # Assume second is _2
          ;;
        *)
          log_info "Warning: Multiple matches for $a, using first two"
          new_fastq1="${matches[0]}"
          new_fastq2="${matches[1]}"
          ;;
      esac
      echo "${b}${DELIMITER}${a}${DELIMITER}${t}${DELIMITER}${new_fastq1}${DELIMITER}${new_fastq2}${DELIMITER}${rest}" >> "$temp_file"
    done < "$INPUT_FILE"

    [[ $missing_count -gt 0 ]] && { log_info "Error: $missing_count samples missing FASTQ files"; rm -f "$temp_file"; exit 1; }

    update_input_csv "$(cat "$temp_file")"
    rm -f "$temp_file"
    log_info "Updated input with FASTQ paths"
  fi
}

########################################
# Initialize/Check Checkpoint Columns
########################################
initialize_checkpoints() {
  log_info "Ensuring all checkpoint columns exist"
  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")

  if [[ "$header" != *",Fastp,Dada2,Motus,Metaphlan,Completed"* ]]; then
    header="${header},Fastp,Dada2,Motus,Metaphlan,Completed"
  fi
  echo "$header" > "$temp_file"

  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c rest || [[ -n "$b" ]]; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    [[ -z "$fp" ]] && fp=0
    [[ -z "$d2" ]] && d2=0
    [[ -z "$mt" ]] && mt=0
    [[ -z "$mp" ]] && mp=0
    [[ -z "$c" ]] && c=0

    echo "${b}${DELIMITER}${a}${DELIMITER}${st}${DELIMITER}${f1}${DELIMITER}${f2}${DELIMITER}${fp}${DELIMITER}${d2}${DELIMITER}${mt}${DELIMITER}${mp}${DELIMITER}${c}" >> "$temp_file"
  done < <(tail -n +2 "$INPUT_FILE")

  update_input_csv "$(cat "$temp_file")"
  rm -f "$temp_file"
  log_info "Checkpoint columns initialized"
}

update_checkpoint() {
  local accession="$1"
  local step="$2"
  local value="$3"
  local field_index=0
  case "$step" in
    Fastp)     field_index=6  ;;
    Dada2)     field_index=7  ;;
    Motus)     field_index=8  ;;
    Metaphlan) field_index=9  ;;
    Completed) field_index=10 ;;
  esac

  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")
  echo "$header" > "$temp_file"

  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c || [[ -n "$b" ]]; do
    if [[ "$b" == "Bioproject" || -z "$b" ]]; then
      continue
    fi
    if [[ "$a" == "$accession" ]]; then
      fields=("$b" "$a" "$st" "$f1" "$f2" "$fp" "$d2" "$mt" "$mp" "$c")
      fields[$((field_index-1))]="$value"
      echo "${fields[0]}${DELIMITER}${fields[1]}${DELIMITER}${fields[2]}${DELIMITER}${fields[3]}${DELIMITER}${fields[4]}${DELIMITER}${fields[5]}${DELIMITER}${fields[6]}${DELIMITER}${fields[7]}${DELIMITER}${fields[8]}${DELIMITER}${fields[9]}" >> "$temp_file"
    else
      echo "$b${DELIMITER}$a${DELIMITER}$st${DELIMITER}$f1${DELIMITER}$f2${DELIMITER}$fp${DELIMITER}$d2${DELIMITER}$mt${DELIMITER}$mp${DELIMITER}$c" >> "$temp_file"
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  update_input_csv "$(cat "$temp_file")"
  rm -f "$temp_file"
  log_debug "Updated checkpoint for $accession: $step=$value"
}

########################################
# Convert Metaphlan profile to counts
########################################
convert_metaphlan_to_counts() {
  local metaphlan_log="$1"
  local metaphlan_profile="$2"
  local output_file="$3"

  if [[ ! -f "$metaphlan_profile" ]]; then
    echo "Missing $metaphlan_profile"
    return 1
  fi

  local mapped_reads
  mapped_reads=$(grep "Total number of reads mapped" "$metaphlan_log" 2>/dev/null | awk '{print $6}' | sed 's/(.*//')

  awk -v mapped="$mapped_reads" '
    BEGIN { FS="\t"; OFS="\t"; print "#clade_name\trelative_abundance\tread_count" }
    /^#/ { next }
    { count = ($2 * mapped / 100); print $1, $2, count }
  ' "$metaphlan_profile" > "$output_file"

  log_debug "Converted $metaphlan_profile to counts in $output_file"
}

########################################
# Merge Profiles
########################################
merge_profiles() {
  local bioproject="$1"
  local tool="$2"
  local output_dir="$bioproject"
  local merged_file="${output_dir}/${bioproject}_${tool}_merged.txt"

  declare -A expected_accessions
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c || [[ -n "$b" ]]; do
    if [[ "$b" == "$bioproject" ]]; then
      expected_accessions["$a"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  local all_complete=true
  for accession in "${!expected_accessions[@]}"; do
    local done_val
    done_val=$(awk -F"$DELIMITER" -v acc="$accession" '$2 == acc {print $10}' "$INPUT_FILE")
    if [[ "$done_val" != "1" ]]; then
      all_complete=false
      break
    fi
  done

  if [[ "$all_complete" == "true" ]]; then
    local profile_files=()
    case "$tool" in
      metaphlan) profile_files=("$output_dir"/*_metaphlan4_counts.txt) ;;
      motus)     profile_files=("$output_dir"/*_motus.txt) ;;
      dada2)     profile_files=("$output_dir"/seqtab_*.rds) ;;
    esac

    local valid_files=()
    for f in "${profile_files[@]}"; do
      [[ -f "$f" ]] && valid_files+=("$f")
    done

    if [[ "$tool" == "metaphlan" && ${#valid_files[@]} -gt 1 ]]; then
      run_command "micromamba run -n \"$MPA_ENV_NAME\" merge_metaphlan_tables.py ${valid_files[*]} > \"$merged_file\"" \
        "[metaphlan merge] for $bioproject" "${valid_files[*]}" "$merged_file"

    elif [[ "$tool" == "metaphlan" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"

    elif [[ "$tool" == "motus" && ${#valid_files[@]} -gt 1 ]]; then
      local input_list
      input_list=$(printf "%s," "${valid_files[@]}" | sed 's/,$//')
      run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus merge -i \"$input_list\" -o \"$merged_file\"" \
        "[motus merge] for $bioproject" "$input_list" "$merged_file"

    elif [[ "$tool" == "motus" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"

    elif [[ "$tool" == "dada2" && ${#valid_files[@]} -gt 0 ]]; then
      local merge_script="$TMP_BASE/merge_dada2.R"
      if [[ ! -f "$merge_script" ]]; then
        log_info "Error: merge_dada2.R not found at $merge_script for DADA2 merging"
        return 1
      fi
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript \"$merge_script\" $bioproject ${valid_files[*]}" \
        "[dada2 merge] for $bioproject" "${valid_files[*]}" ""
    fi
  fi
}


########################################
# process_sample function
########################################
process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SEQ_TYPE="$3"
  local FASTQ1="$4"
  local FASTQ2="$5"
  local FASTP="$6"
  local DADA2="$7"
  local MOTUS="$8"
  local METAPHLAN="$9"
  local COMPLETED="${10}"

  log_info "Starting process_sample for $RUN_ACCESSION"

  if [[ "$COMPLETED" == "1" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already completed"
    return
  fi

  local OUTPUT_DIR="${BIOPROJECT}"
  log_info "Creating output directory for $RUN_ACCESSION: $OUTPUT_DIR"
  mkdir -p "$OUTPUT_DIR"
  local TMP_DIR
  TMP_DIR=$(mktemp -d -p "$TMP_BASE" "process_${SLURM_JOB_ID}_${RUN_ACCESSION}_XXXXXX")
  TMP_DIRS+=("$TMP_DIR")
  log_info "Created temporary directory for $RUN_ACCESSION: $TMP_DIR"

  local QC1="${TMP_DIR}/qc_${RUN_ACCESSION}_1.fastq.gz"
  local QC2="${TMP_DIR}/qc_${RUN_ACCESSION}_2.fastq.gz"
  local QC="${TMP_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
  local METAPHLAN_LOG="${TMP_DIR}/${RUN_ACCESSION}_metaphlan_log.txt"
  local METAPHLAN_BOWTIE="${TMP_DIR}/${RUN_ACCESSION}_meta.bowtie2out.txt"
  local METAPHLAN_PROFILE="${TMP_DIR}/${RUN_ACCESSION}_metaphlan4.txt"
  local METAPHLAN_COUNTS="${OUTPUT_DIR}/${RUN_ACCESSION}_metaphlan4_counts.txt"
  local MOTUS_PROFILE="${OUTPUT_DIR}/${RUN_ACCESSION}_motus.txt"
  local DADA2_SEQTAB="${OUTPUT_DIR}/seqtab_${RUN_ACCESSION}.rds"

  log_info "Validating FASTQ files for $RUN_ACCESSION"
  if [[ ! -f "$FASTQ1" ]]; then
    append_with_lock "$RUN_ACCESSION:Missing FASTQ $FASTQ1" "$FAILED_FILE" "$FAILED_LOCK"
    return 1
  else
    local reason1
    reason1=$(validate_fastq "$FASTQ1") || {
      append_with_lock "$RUN_ACCESSION:$reason1" "$FAILED_FILE" "$FAILED_LOCK"
      return 1
    }
  fi

  if [[ -n "$FASTQ2" && -f "$FASTQ2" ]]; then
    local reason2
    reason2=$(validate_fastq "$FASTQ2") || {
      append_with_lock "$RUN_ACCESSION:$reason2" "$FAILED_FILE" "$FAILED_LOCK"
      return 1
    }
  fi
  log_info "Finished validating FASTQ for $RUN_ACCESSION"

  local FASTQ_TO_USE_1="$FASTQ1"
  local FASTQ_TO_USE_2="$FASTQ2"
  if [[ "$QC_ENABLED" == "true" && "$FASTP" != "1" ]]; then
    log_info "Starting fastp for $RUN_ACCESSION"
    if [[ -n "$FASTQ2" ]]; then
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$FASTQ1\" -I \"$FASTQ2\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER -Q -A -L" \
        "[fastp] Process for $RUN_ACCESSION" "$FASTQ1,$FASTQ2" "$QC1,$QC2"
      if [[ $? -eq 0 && -f "$QC1" && -f "$QC2" ]]; then
        local size1 size2
        size1=$(stat -c %s "$QC1" 2>/dev/null || wc -c < "$QC1")
        size2=$(stat -c %s "$QC2" 2>/dev/null || wc -c < "$QC2")
        if [[ $size1 -gt 100 && $size2 -gt 100 ]]; then
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
          FASTQ_TO_USE_1="$QC1"
          FASTQ_TO_USE_2="$QC2"
        else
          log_info "Error: fastp produced empty or too small outputs for $RUN_ACCESSION ($size1, $size2 bytes)"
          append_with_lock "$RUN_ACCESSION:fastp produced empty outputs ($size1, $size2 bytes)" "$FAILED_FILE" "$FAILED_LOCK"
          return 1
        fi
      else
        log_info "Error: fastp failed or missing outputs for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:fastp failed or missing $QC1/$QC2" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
    else
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$FASTQ1\" -o \"$QC\" -w $THREADS_PER_WORKER -Q -A -L" \
        "[fastp] Process for $RUN_ACCESSION" "$FASTQ1" "$QC"
      if [[ $? -eq 0 && -f "$QC" ]]; then
        local size
        size=$(stat -c %s "$QC" 2>/dev/null || wc -c < "$QC")
        if [[ $size -gt 100 ]]; then
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
          FASTQ_TO_USE_1="$QC"
        else
          log_info "Error: fastp produced empty or too small output for $RUN_ACCESSION ($size bytes)"
          append_with_lock "$RUN_ACCESSION:fastp produced empty output ($size bytes)" "$FAILED_FILE" "$FAILED_LOCK"
          return 1
        fi
      else
        log_info "Error: fastp failed or missing output for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:fastp failed or missing $QC" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
    fi
    log_info "Finished fastp for $RUN_ACCESSION"
  else
    log_debug "Skipping fastp for $RUN_ACCESSION"
  fi

  if [[ "$SEQ_TYPE" == "16S" && "$DADA2" != "1" ]]; then
    log_info "Starting DADA2 for $RUN_ACCESSION"
    if [[ -n "$FASTQ_TO_USE_2" ]]; then
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript \"$TMP_BASE/run_dada2_partial.R\" \"$FASTQ_TO_USE_1\" \"$FASTQ_TO_USE_2\")" \
        "[dada2] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1,$FASTQ_TO_USE_2" "$DADA2_SEQTAB"
      if [[ $? -eq 0 && -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
        if [[ -f "$DADA2_SEQTAB" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Dada2" "1"
        else
          log_info "Error: DADA2 output move failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:DADA2 output move failed" "$FAILED_FILE" "$FAILED_LOCK"
          return 1
        fi
      else
        log_info "Error: DADA2 failed or output not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 failed or missing output" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
    else
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript \"$TMP_BASE/run_dada2_partial.R\" \"$FASTQ_TO_USE_1\")" \
        "[dada2] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1" "$DADA2_SEQTAB"
      if [[ $? -eq 0 && -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
        if [[ -f "$DADA2_SEQTAB" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Dada2" "1"
        else
          log_info "Error: DADA2 output move failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:DADA2 output move failed" "$FAILED_FILE" "$FAILED_LOCK"
          return 1
        fi
      else
        log_info "Error: DADA2 failed or output not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 failed or missing output" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
    fi
    log_info "Finished DADA2 for $RUN_ACCESSION"
  else
    log_debug "Skipping DADA2 for $RUN_ACCESSION"
  fi

  if [[ "$SEQ_TYPE" == "meta" ]]; then
    if [[ "$METAPHLAN" != "1" ]]; then
      log_info "Starting MetaPhlAn for $RUN_ACCESSION"
      local FASTQ_INPUT="$FASTQ_TO_USE_1"
      [[ -n "$FASTQ_TO_USE_2" ]] && FASTQ_INPUT="$FASTQ_TO_USE_1,$FASTQ_TO_USE_2"
      run_command "micromamba run -n \"$MPA_ENV_NAME\" metaphlan \"$FASTQ_INPUT\" --input_type fastq --unclassified_estimation --nproc "$THREADS_PER_WORKER" --bowtie2db \"$METAPHLAN_DB\" --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\" 2> \"$METAPHLAN_LOG\"" \
        "[metaphlan] Process for $RUN_ACCESSION" "$FASTQ_INPUT" "$METAPHLAN_PROFILE"
      if [[ $? -eq 0 && -f "$METAPHLAN_PROFILE" ]]; then
        convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS"
        if [[ -f "$METAPHLAN_COUNTS" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Metaphlan" "1"
        else
          log_info "Error: MetaPhlAn counts conversion failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:MetaPhlAn counts conversion failed" "$FAILED_FILE" "$FAILED_LOCK"
          return 1
        fi
      else
        local metaphlan_error="Unknown error"
        [[ -f "$METAPHLAN_LOG" ]] && metaphlan_error=$(tail -n 5 "$METAPHLAN_LOG" | tr '\n' ';')
        log_info "Error: MetaPhlAn failed for $RUN_ACCESSION. Log tail: $metaphlan_error"
        append_with_lock "$RUN_ACCESSION:MetaPhlAn failed. Log tail: $metaphlan_error" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
      log_info "Finished MetaPhlAn for $RUN_ACCESSION"
    else
      log_debug "Skipping MetaPhlAn for $RUN_ACCESSION"
    fi

    if [[ "$MOTUS" != "1" ]]; then
      log_info "Starting mOTUs for $RUN_ACCESSION"
      if [[ -n "$FASTQ_TO_USE_2" ]]; then
        run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -f \"$FASTQ_TO_USE_1\" -r \"$FASTQ_TO_USE_2\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1,$FASTQ_TO_USE_2" "$MOTUS_PROFILE"
      else
        run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -s \"$FASTQ_TO_USE_1\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1" "$MOTUS_PROFILE"
      fi
      if [[ $? -eq 0 && -f "$MOTUS_PROFILE" ]]; then
        update_checkpoint "$RUN_ACCESSION" "Motus" "1"
      else
        log_info "Error: mOTUs failed or missing profile for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:mOTUs failed or missing profile" "$FAILED_FILE" "$FAILED_LOCK"
        return 1
      fi
      log_info "Finished mOTUs for $RUN_ACCESSION"
    else
      log_debug "Skipping mOTUs for $RUN_ACCESSION"
    fi
  fi

  if [[ "$SEQ_TYPE" == "16S" && -f "$DADA2_SEQTAB" ]] || \
     [[ "$SEQ_TYPE" == "meta" && -f "$METAPHLAN_COUNTS" && -f "$MOTUS_PROFILE" ]]; then
    update_checkpoint "$RUN_ACCESSION" "Completed" "1"
    log_info "Completed processing for $RUN_ACCESSION"
  else
    log_info "Error: Not all outputs verified for $RUN_ACCESSION"
    append_with_lock "$RUN_ACCESSION:Incomplete outputs" "$FAILED_FILE" "$FAILED_LOCK"
    return 1
  fi
}

# Export functions and variables for parallel
export -f process_sample log_debug log_info cleanup update_input_csv append_with_lock update_checkpoint run_command run_command_with_output convert_metaphlan_to_counts validate_fastq
export DELIMITER INPUT_FILE FAILED_FILE FAILED_LOCK INPUT_LOCK FASTQ_DIR MPA_ENV_NAME METAPHLAN_DB THREADS_PER_WORKER DADA2_ENV_NAME MOTUS_ENV_NAME MOTUS_TAX_LEVEL TMP_BASE LOG_LEVEL

########################################
# Debug Mode Verification
########################################
verify_sample_output() {
  local accession="$1"
  local seq_type="$2"
  local output_dir="$3"
  local status=0

  if [[ "$seq_type" == "16S" ]]; then
    local dada2_seqtab="${output_dir}/seqtab_${accession}.rds"
    if [[ ! -f "$dada2_seqtab" ]]; then
      log_info "Debug verification failed for $accession: Missing DADA2 output $dada2_seqtab"
      status=1
    else
      log_info "Debug verification passed for $accession: DADA2 output $dada2_seqtab exists"
    fi
  elif [[ "$seq_type" == "meta" ]]; then
    local metaphlan_counts="${output_dir}/${accession}_metaphlan4_counts.txt"
    local motus_profile="${output_dir}/${accession}_motus.txt"
    if [[ ! -f "$metaphlan_counts" ]]; then
      log_info "Debug verification failed for $accession: Missing MetaPhlAn output $metaphlan_counts"
      status=1
    else
      log_info "Debug verification passed for $accession: MetaPhlAn output $metaphlan_counts exists"
    fi
    if [[ ! -f "$motus_profile" ]]; then
      log_info "Debug verification failed for $accession: Missing mOTUs output $motus_profile"
      status=1
    else
      log_info "Debug verification passed for $accession: mOTUs output $motus_profile exists"
    fi
  fi
  return $status
}

########################################
# Final Validation & Merging
########################################
final_validation_and_merge() {
  log_info "Starting final validation and merging"

  declare -A BIOPROJECTS
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c || [[ -n "$b" ]]; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    BIOPROJECTS["$b"]=1
  done < <(tail -n +2 "$INPUT_FILE")

  for bioproject in "${!BIOPROJECTS[@]}"; do
    merge_profiles "$bioproject" "metaphlan"
    merge_profiles "$bioproject" "motus"
    merge_profiles "$bioproject" "dada2"
  done

  log_info "Final validation and merging complete."
}

########################################
# Debug Mode Parallel Processing
########################################
debug_mode_run() {
  log_info "DEBUG mode: Processing one 16S and one meta sample in parallel"

  local debug_input
  debug_input=$(mktemp)
  local processed_16s=0
  local processed_meta=0

  head -n 1 "$INPUT_FILE" > "$debug_input"

  dos2unix 2>/dev/null "$INPUT_FILE" || true
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c || [[ -n "$b" ]]; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    st="$(echo "$st" | tr -d '\r')"

    if [[ "$st" == "16S" && $processed_16s -eq 0 ]]; then
      echo "$b${DELIMITER}$a${DELIMITER}$st${DELIMITER}$f1${DELIMITER}$f2${DELIMITER}$fp${DELIMITER}$d2${DELIMITER}$mt${DELIMITER}$mp${DELIMITER}$c" >> "$debug_input"
      ((processed_16s++))
    elif [[ "$st" == "meta" && $processed_meta -eq 0 ]]; then
      echo "$b${DELIMITER}$a${DELIMITER}$st${DELIMITER}$f1${DELIMITER}$f2${DELIMITER}$fp${DELIMITER}$d2${DELIMITER}$mt${DELIMITER}$mp${DELIMITER}$c" >> "$debug_input"
      ((processed_meta++))
    fi

    if [[ $processed_16s -eq 1 && $processed_meta -eq 1 ]]; then
      break
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  if [[ $processed_16s -eq 0 || $processed_meta -eq 0 ]]; then
    log_info "Error: Could not find both a 16S and a meta sample for debug mode"
    rm -f "$debug_input"
    exit 1
  fi

  # NEW: Redirect parallel errors to a log file for debugging
  local parallel_error_log
  parallel_error_log=$(mktemp)
  tail -n +2 "$debug_input" | \
    parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" --halt now,fail=1 \
      --env process_sample \
      --env log_debug --env log_info --env append_with_lock --env update_checkpoint \
      --env run_command --env run_command_with_output --env convert_metaphlan_to_counts \
      --env validate_fastq --env update_input_csv \
      --env DELIMITER --env INPUT_FILE --env FAILED_FILE --env FAILED_LOCK \
      --env INPUT_LOCK \
      --env FASTQ_DIR --env MPA_ENV_NAME --env METAPHLAN_DB \
      --env THREADS_PER_WORKER --env DADA2_ENV_NAME --env MOTUS_ENV_NAME \
      --env MOTUS_TAX_LEVEL --env TMP_BASE --env LOG_LEVEL \
      "bash -c 'process_sample {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'" 2> "$parallel_error_log"
  if [[ $? -ne 0 ]]; then
    echo "Error: Parallel processing failed in debug mode: $(cat "$parallel_error_log")"
    rm -f "$debug_input" "$parallel_error_log"
    exit 1
  fi
  rm -f "$parallel_error_log"

  processed_16s=0
  processed_meta=0
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c || [[ -n "$b" ]]; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    if [[ "$st" == "16S" && $processed_16s -eq 0 ]]; then
      if verify_sample_output "$a" "$st" "$b"; then
        log_info "Debug mode: 16S sample $a processed successfully"
      else
        log_info "Debug mode: 16S sample $a failed verification"
        rm -f "$debug_input"
        exit 1
      fi
      ((processed_16s++))
    elif [[ "$st" == "meta" && $processed_meta -eq 0 ]]; then
      if verify_sample_output "$a" "$st" "$b"; then
        log_info "Debug mode: Meta sample $a processed successfully"
      else
        log_info "Debug mode: Meta sample $a failed verification"
        rm -f "$debug_input"
        exit 1
      fi
      ((processed_meta++))
    fi
  done < <(tail -n +2 "$debug_input")

  rm -f "$debug_input"
  log_info "DEBUG mode: Done processing 1 16S and 1 meta sample"
}

########################################
# Main Workflow
########################################
update_input_with_fastq_paths
initialize_checkpoints

log_info "Starting processing with $NUM_WORKERS workers"

if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
  debug_mode_run
else
  # NEW: Redirect parallel errors to a log file for debugging
  local parallel_error_log
  parallel_error_log=$(mktemp)
  tail -n +2 "$INPUT_FILE" | \
    parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" --halt now,fail=1 \
      --env process_sample \
      --env log_debug --env log_info --env append_with_lock --env update_checkpoint \
      --env run_command --env run_command_with_output --env convert_metaphlan_to_counts \
      --env validate_fastq --env update_input_csv \
      --env DELIMITER --env INPUT_FILE --env FAILED_FILE --env FAILED_LOCK \
      --env INPUT_LOCK \
      --env FASTQ_DIR --env MPA_ENV_NAME --env METAPHLAN_DB \
      --env THREADS_PER_WORKER --env DADA2_ENV_NAME --env MOTUS_ENV_NAME \
      --env MOTUS_TAX_LEVEL --env TMP_BASE --env LOG_LEVEL \
      "bash -c 'process_sample {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'" 2> "$parallel_error_log"
  if [[ $? -ne 0 ]]; then
    echo "Error: Parallel processing failed: $(cat "$parallel_error_log")"
    rm -f "$parallel_error_log"
    exit 1
  fi
  rm -f "$parallel_error_log"
fi

final_validation_and_merge

log_info "All processing, validation, and merging complete."
echo "All processing, validation, and merging complete."
########################################
# Initialize/Check Checkpoint Columns
########################################
initialize_checkpoints() {
  log_info "Ensuring all checkpoint columns exist"
  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")

  if [[ "$header" != *",Fastp,Dada2,Motus,Metaphlan,Completed"* ]]; then
    header="${header},Fastp,Dada2,Motus,Metaphlan,Completed"
  fi
  echo "$header" > "$temp_file"

  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c rest; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    [[ -z "$fp" ]] && fp=0
    [[ -z "$d2" ]] && d2=0
    [[ -z "$mt" ]] && mt=0
    [[ -z "$mp" ]] && mp=0
    [[ -z "$c" ]] && c=0

    echo "${b}${DELIMITER}${a}${DELIMITER}${st}${DELIMITER}${f1}${DELIMITER}${f2}${DELIMITER}${fp}${DELIMITER}${d2}${DELIMITER}${mt}${DELIMITER}${mp}${DELIMITER}${c}" >> "$temp_file"
  done < <(tail -n +2 "$INPUT_FILE")

  update_input_csv "$(cat "$temp_file")"
  rm -f "$temp_file"
  log_info "Checkpoint columns initialized"
}

update_checkpoint() {
  local accession="$1"
  local step="$2"
  local value="$3"
  local field_index=0
  case "$step" in
    Fastp)     field_index=6  ;;
    Dada2)     field_index=7  ;;
    Motus)     field_index=8  ;;
    Metaphlan) field_index=9  ;;
    Completed) field_index=10 ;;
  esac

  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")
  echo "$header" > "$temp_file"

  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c; do
    if [[ "$b" == "Bioproject" || -z "$b" ]]; then
      continue
    fi
    if [[ "$a" == "$accession" ]]; then
      fields=("$b" "$a" "$st" "$f1" "$f2" "$fp" "$d2" "$mt" "$mp" "$c")
      fields[$((field_index-1))]="$value"
      echo "${fields[0]}${DELIMITER}${fields[1]}${DELIMITER}${fields[2]}${DELIMITER}${fields[3]}${DELIMITER}${fields[4]}${DELIMITER}${fields[5]}${DELIMITER}${fields[6]}${DELIMITER}${fields[7]}${DELIMITER}${fields[8]}${DELIMITER}${fields[9]}" >> "$temp_file"
    else
      echo "$b${DELIMITER}$a${DELIMITER}$st${DELIMITER}$f1${DELIMITER}$f2${DELIMITER}$fp${DELIMITER}$d2${DELIMITER}$mt${DELIMITER}$mp${DELIMITER}$c" >> "$temp_file"
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  update_input_csv "$(cat "$temp_file")"
  rm -f "$temp_file"
  log_debug "Updated checkpoint for $accession: $step=$value"
}

########################################
# Convert Metaphlan profile to counts
########################################
convert_metaphlan_to_counts() {
  local metaphlan_log="$1"
  local metaphlan_profile="$2"
  local output_file="$3"

  if [[ ! -f "$metaphlan_profile" ]]; then
    echo "Missing $metaphlan_profile"
    return 1
  fi

  local mapped_reads
  mapped_reads=$(grep "Total number of reads mapped" "$metaphlan_log" 2>/dev/null | awk '{print $6}' | sed 's/(.*//')
  if [[ -z "$mapped_reads" ]]; then
    mapped_reads=100000
  fi

  awk -v mapped="$mapped_reads" '
    BEGIN { FS="\t"; OFS="\t"; print "#clade_name\trelative_abundance\tread_count" }
    /^#/ { next }
    { count = ($2 * mapped / 100); print $1, $2, count }
  ' "$metaphlan_profile" > "$output_file"

  log_debug "Converted $metaphlan_profile to counts in $output_file"
}

########################################
# Merge Profiles
########################################
merge_profiles() {
  local bioproject="$1"
  local tool="$2"
  local output_dir="$bioproject"
  local merged_file="${output_dir}/${bioproject}_${tool}_merged.txt"

  declare -A expected_accessions
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c; do
    if [[ "$b" == "$bioproject" ]]; then
      expected_accessions["$a"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  local all_complete=true
  for accession in "${!expected_accessions[@]}"; do
    local done_val
    done_val=$(awk -F"$DELIMITER" -v acc="$accession" '$2 == acc {print $10}' "$INPUT_FILE")
    if [[ "$done_val" != "1" ]]; then
      all_complete=false
      break
    fi
  done

  if [[ "$all_complete" == "true" ]]; then
    local profile_files=()
    case "$tool" in
      metaphlan)   profile_files=("$output_dir"/*_metaphlan4_counts.txt) ;;
      motus)       profile_files=("$output_dir"/*_motus.txt) ;;
      dada2)       profile_files=("$output_dir"/seqtab_*.rds) ;;
    esac

    local valid_files=()
    for f in "${profile_files[@]}"; do
      [[ -f "$f" ]] && valid_files+=("$f")
    done

    if [[ "$tool" == "metaphlan" && ${#valid_files[@]} -gt 1 ]]; then
      run_command "micromamba run -n \"$MPA_ENV_NAME\" merge_metaphlan_tables.py ${valid_files[*]} > \"$merged_file\"" \
        "[metaphlan merge] for $bioproject" "${valid_files[*]}" "$merged_file"

    elif [[ "$tool" == "metaphlan" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"

    elif [[ "$tool" == "motus" && ${#valid_files[@]} -gt 1 ]]; then
      local input_list
      input_list=$(printf "%s," "${valid_files[@]}" | sed 's/,$//')
      run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus merge -i \"$input_list\" -o \"$merged_file\"" \
        "[motus merge] for $bioproject" "$input_list" "$merged_file"

    elif [[ "$tool" == "motus" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"

    elif [[ "$tool" == "dada2" && ${#valid_files[@]} -gt 0 ]]; then
      if [[ ! -f "merge_dada2.R" ]]; then
        log_info "Error: merge_dada2.R not found for DADA2 merging"
        return 1
      fi
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript merge_dada2.R $bioproject ${valid_files[*]}" \
        "[dada2 merge] for $bioproject" "${valid_files[*]}" ""
    fi
  fi
}

########################################
# Global Array for Completed Samples
########################################
declare -A COMPLETED_SAMPLES
load_completed() {
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c; do
    if [[ "$c" == "1" ]]; then
      COMPLETED_SAMPLES["$a"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")
  log_debug "Loaded ${#COMPLETED_SAMPLES[@]} completed samples"
}

########################################
# process_sample function
########################################
process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SAMPLE_TYPE="$3"
  local INPUT_FASTQ="$4"
  local PAIRED_FASTQ="$5"

  log_debug "Processing $RUN_ACCESSION ($SAMPLE_TYPE)"

  if [[ -n "${COMPLETED_SAMPLES[$RUN_ACCESSION]}" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already completed"
    return
  fi

  local OUTPUT_DIR="${BIOPROJECT}"
  mkdir -p "$OUTPUT_DIR"
  local TMP_DIR
  TMP_DIR=$(mktemp -d -p "$TMP_BASE" "process_${SLURM_JOB_ID}_${RUN_ACCESSION}_XXXXXX")
  TMP_DIRS+=("$TMP_DIR")

  local QC1="${TMP_DIR}/qc_${RUN_ACCESSION}_1.fastq.gz"
  local QC2="${TMP_DIR}/qc_${RUN_ACCESSION}_2.fastq.gz"
  local QC="${TMP_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
  local METAPHLAN_LOG="${TMP_DIR}/${RUN_ACCESSION}_metaphlan_log.txt"
  local METAPHLAN_BOWTIE="${TMP_DIR}/${RUN_ACCESSION}_meta.bowtie2out.txt"
  local METAPHLAN_PROFILE="${TMP_DIR}/${RUN_ACCESSION}_metaphlan4.txt"
  local METAPHLAN_COUNTS="${OUTPUT_DIR}/${RUN_ACCESSION}_metaphlan4_counts.txt"
  local MOTUS_PROFILE="${OUTPUT_DIR}/${RUN_ACCESSION}_motus.txt"
  local DADA2_SEQTAB="${OUTPUT_DIR}/seqtab_${RUN_ACCESSION}.rds"

  # Validate input fastq(s)
  if [[ ! -f "$INPUT_FASTQ" ]]; then
    append_with_lock "$RUN_ACCESSION:Missing FASTQ $INPUT_FASTQ" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  else
    local reason1
    reason1=$(validate_fastq "$INPUT_FASTQ") || {
      append_with_lock "$RUN_ACCESSION:$reason1" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    }
  fi

  if [[ -n "$PAIRED_FASTQ" && -f "$PAIRED_FASTQ" ]]; then
    local reason2
    reason2=$(validate_fastq "$PAIRED_FASTQ") || {
      append_with_lock "$RUN_ACCESSION:$reason2" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    }
  fi

  # Decide whether to run fastp
  local FASTQ_TO_USE_1="$INPUT_FASTQ"
  local FASTQ_TO_USE_2="$PAIRED_FASTQ"
  if [[ "$QC_ENABLED" == "true" ]]; then
    log_info "Running fastp QC for $RUN_ACCESSION"
    local fastp_status
    fastp_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $6}' "$INPUT_FILE")
    if [[ "$fastp_status" != "1" ]]; then
      if [[ -n "$PAIRED_FASTQ" ]]; then
        run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER -Q -A -L" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2"
        if [[ $? -eq 0 && -f "$QC1" && -f "$QC2" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
          FASTQ_TO_USE_1="$QC1"
          FASTQ_TO_USE_2="$QC2"
        else
          log_info "Error: fastp failed or missing outputs for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:fastp failed or missing $QC1/$QC2" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      else
        run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER -Q -A -L" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC"
        if [[ $? -eq 0 && -f "$QC" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
          FASTQ_TO_USE_1="$QC"
        else
          log_info "Error: fastp failed or missing output for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:fastp failed or missing $QC" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
    else
      log_debug "Skipping fastp for $RUN_ACCESSION: checkpoint set"
      FASTQ_TO_USE_1="$QC1"
      FASTQ_TO_USE_2="$QC2"
    fi
  else
    log_info "Skipping fastp for $RUN_ACCESSION"
  fi

  ################################
  # 16S Workflow
  ################################
  if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    if [[ -n "$FASTQ_TO_USE_2" ]]; then
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$FASTQ_TO_USE_1\" \"$FASTQ_TO_USE_2\")" \
        "[dada2] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1,$FASTQ_TO_USE_2" "$DADA2_SEQTAB"
      if [[ $? -eq 0 && -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
        if [[ -f "$DADA2_SEQTAB" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Dada2" "1"
        else
          log_info "Error: DADA2 output move failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:DADA2 output move failed" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      else
        log_info "Error: DADA2 failed or output not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 failed or missing output" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    else
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$FASTQ_TO_USE_1\")" \
        "[dada2] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1" "$DADA2_SEQTAB"
      if [[ $? -eq 0 && -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
        if [[ -f "$DADA2_SEQTAB" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Dada2" "1"
        else
          log_info "Error: DADA2 output move failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:DADA2 output move failed" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      else
        log_info "Error: DADA2 failed or output not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 failed or missing output" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    fi

  ################################
  # Metagenomic Workflow
  ################################
  elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    local metaphlan_status
    metaphlan_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $9}' "$INPUT_FILE")
    if [[ "$metaphlan_status" != "1" ]]; then
      local FASTQ_INPUT="$FASTQ_TO_USE_1"
      [[ -n "$FASTQ_TO_USE_2" ]] && FASTQ_INPUT="$FASTQ_TO_USE_1,$FASTQ_TO_USE_2"
      run_command "micromamba run -n \"$MPA_ENV_NAME\" metaphlan \"$FASTQ_INPUT\" --input_type fastq --unclassified_estimation --nproc $THREADS_PER_WORKER --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\" 2> \"$METAPHLAN_LOG\"" \
        "[metaphlan] Process for $RUN_ACCESSION" "$FASTQ_INPUT" "$METAPHLAN_PROFILE"
      if [[ $? -eq 0 && -f "$METAPHLAN_PROFILE" ]]; then
        convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS"
        if [[ -f "$METAPHLAN_COUNTS" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Metaphlan" "1"
        else
          log_info "Error: MetaPhlAn counts conversion failed for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:MetaPhlAn counts conversion failed" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      else
        log_info "Error: MetaPhlAn failed or missing profile for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:MetaPhlAn failed or missing profile" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    else
      log_debug "Skipping metaphlan for $RUN_ACCESSION: checkpoint set"
    fi

    local motus_status
    motus_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $8}' "$INPUT_FILE")
    if [[ "$motus_status" != "1" ]]; then
      if [[ -n "$FASTQ_TO_USE_2" ]]; then
        run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -f \"$FASTQ_TO_USE_1\" -r \"$FASTQ_TO_USE_2\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1,$FASTQ_TO_USE_2" "$MOTUS_PROFILE"
        if [[ $? -eq 0 && -f "$MOTUS_PROFILE" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Motus" "1"
        else
          log_info "Error: mOTUs failed or missing profile for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:mOTUs failed or missing profile" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      else
        run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -s \"$FASTQ_TO_USE_1\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$FASTQ_TO_USE_1" "$MOTUS_PROFILE"
        if [[ $? -eq 0 && -f "$MOTUS_PROFILE" ]]; then
          update_checkpoint "$RUN_ACCESSION" "Motus" "1"
        else
          log_info "Error: mOTUs failed or missing profile for $RUN_ACCESSION"
          append_with_lock "$RUN_ACCESSION:mOTUs failed or missing profile" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
    else
      log_debug "Skipping motus for $RUN_ACCESSION: checkpoint set"
    fi
  else
    local reason="Invalid sample type: $SAMPLE_TYPE"
    append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  # Only mark as completed if all relevant outputs are verified
  if [[ "$SAMPLE_TYPE" == "16S" && -f "$DADA2_SEQTAB" ]] || \
     [[ "$SAMPLE_TYPE" == "meta" && -f "$METAPHLAN_COUNTS" && -f "$MOTUS_PROFILE" ]]; then
    update_checkpoint "$RUN_ACCESSION" "Completed" "1"
    append_with_lock "${RUN_ACCESSION}:COMPLETE" "$COMPLETED_FILE" "$COMPLETED_LOCK"
    log_info "Finished $RUN_ACCESSION"
  else
    log_info "Error: Not all outputs verified for $RUN_ACCESSION"
    append_with_lock "$RUN_ACCESSION:Incomplete outputs" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  rm -rf "$TMP_DIR"
}

########################################
# Debug Mode Verification
########################################
verify_sample_output() {
  local accession="$1"
  local sample_type="$2"
  local output_dir="$3"
  local status=0

  if [[ "$sample_type" == "16S" ]]; then
    local dada2_seqtab="${output_dir}/seqtab_${accession}.rds"
    if [[ ! -f "$dada2_seqtab" ]]; then
      log_info "Debug verification failed for $accession: Missing DADA2 output $dada2_seqtab"
      status=1
    else
      log_info "Debug verification passed for $accession: DADA2 output $dada2_seqtab exists"
    fi
  elif [[ "$sample_type" == "meta" ]]; then
    local metaphlan_counts="${output_dir}/${accession}_metaphlan4_counts.txt"
    local motus_profile="${output_dir}/${accession}_motus.txt"
    if [[ ! -f "$metaphlan_counts" ]]; then
      log_info "Debug verification failed for $accession: Missing MetaPhlAn output $metaphlan_counts"
      status=1
    else
      log_info "Debug verification passed for $accession: MetaPhlAn output $metaphlan_counts exists"
    fi
    if [[ ! -f "$motus_profile" ]]; then
      log_info "Debug verification failed for $accession: Missing mOTUs output $motus_profile"
      status=1
    else
      log_info "Debug verification passed for $accession: mOTUs output $motus_profile exists"
    fi
  fi
  return $status
}

########################################
# Final Validation & Merging
########################################
final_validation_and_merge() {
  log_info "Starting final validation and merging"

  declare -A BIOPROJECTS
  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    BIOPROJECTS["$b"]=1
  done < <(tail -n +2 "$INPUT_FILE")

  for bioproject in "${!BIOPROJECTS[@]}"; do
    merge_profiles "$bioproject" "metaphlan"
    merge_profiles "$bioproject" "motus"
    merge_profiles "$bioproject" "dada2"
  done

  log_info "Final validation and merging complete."
  log_info "All processing complete."
}

########################################
# Main Workflow
########################################
update_input_with_fastq_paths
initialize_checkpoints
load_completed

debug_mode_run() {
  local processed_16s=0
  local processed_meta=0

  dos2unix 2>/dev/null "$INPUT_FILE" || true

  while IFS="$DELIMITER" read -r b a st f1 f2 fp d2 mt mp c; do
    if [[ -z "$b" || "$b" == "Bioproject" ]]; then
      continue
    fi
    st="$(echo "$st" | tr -d '\r')"

    if [[ "$st" == "16S" && $processed_16s -eq 0 ]]; then
      log_info "Debug mode: Processing 16S sample $a"
      process_sample "$b" "$a" "$st" "$f1" "$f2"
      if verify_sample_output "$a" "$st" "$b"; then
        log_info "Debug mode: 16S sample $a processed successfully"
      else
        log_info "Debug mode: 16S sample $a failed verification"
        exit 1
      fi
      ((processed_16s++))
    elif [[ "$st" == "meta" && $processed_meta -eq 0 ]]; then
      log_info "Debug mode: Processing meta sample $a"
      process_sample "$b" "$a" "$st" "$f1" "$f2"
      if verify_sample_output "$a" "$st" "$b"; then
        log_info "Debug mode: Meta sample $a processed successfully"
      else
        log_info "Debug mode: Meta sample $a failed verification"
        exit 1
      fi
      ((processed_meta++))
    fi

    if [[ $processed_16s -eq 1 && $processed_meta -eq 1 ]]; then
      log_info "DEBUG mode: Done processing 1 16S and 1 meta."
      break
    fi
  done < <(tail -n +2 "$INPUT_FILE")
}

log_info "Starting parallel processing with $NUM_WORKERS workers"

if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
  log_info "DEBUG mode: Process only ONE 16S and ONE meta sample (sequentially)."
  debug_mode_run
else
  tail -n +2 "$INPUT_FILE" | \
    parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" --halt now,fail=1 \
      "$0" --process-sample {1} {2} {3} {4} {5}
fi

final_validation_and_merge

log_info "All processing, validation, and merging complete."
echo "All processing, validation, and merging complete."

# Handle inline process-sample command
if [[ "$1" == "--process-sample" ]]; then
  shift
  process_sample "$@"
  exit $?
fi
