#!/bin/bash

#SBATCH --export=ALL
#SBATCH --job-name=main_concurrent
#SBATCH --output=slurm_%A.out
#SBATCH --error=slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

########################################
# Default Config
########################################
ITEMS_FILE="run_accessions.txt"
COMPLETED_FILE="completed_accessions.txt"
MAX_JOBS=30  # Explicitly set maximum concurrent jobs
JOB_NAME_PREFIX="concurrent_"
JOB_COMMAND="python3 /path/to/download_accession.py"  # Update this path
SLURM_TIME="02:00:00"
SLURM_MEM="8G"
CPUS="4"
THIS_JOB_ID="${SLURM_JOB_ID:-$$}"
ORIGINAL_COMMAND="sbatch $0 $@"
START_TIME="$(date +%s)"
MAX_RUNTIME=$(( 47 * 3600 - 60 ))
LAST_RESET_TIME=0
RESET_INTERVAL=600  
NCBI_API_KEY=""
EMAIL=""
USER_EXPORTS=""
DYNAMIC_RESOURCES=0  
WORKDIR="$(pwd)"
LOG_DIR="${WORKDIR}/logs"
JOBS_DIR="${WORKDIR}/jobs"
TOKEN_FILE="${LOG_DIR}/.job_tokens"
TOKEN_LOCK_FILE="${LOG_DIR}/.job_tokens.lock"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
CHECKPOINT_LOCK_FILE="${LOG_DIR}/checkpoint.lock"
DEBUG_LOCK_FILE="${LOG_DIR}/debug.lock"
PERIODIC_RESET_LOG="${LOG_DIR}/periodic_reset.log"
CLEANUP_COUNTER=0

########################################
# Usage function
########################################
usage() {
  echo "Usage: $0 [options]"
  echo "  -i <file>        Input file of items/tasks (default: $ITEMS_FILE)"
  echo "  -o <file>        Completed items file (default: $COMPLETED_FILE)"
  echo "  -m <int>         Max concurrency (default: $MAX_JOBS)"
  echo "  -p <prefix>      Slurm job-name prefix (default: $JOB_NAME_PREFIX)"
  echo "  -C <cmd>         Command to run per item (default: \"$JOB_COMMAND\")"
  echo "  -T <time>        Slurm time limit per job (default: $SLURM_TIME)"
  echo "  -M <mem>         Slurm memory per job (default: $SLURM_MEM)"
  echo "  --email <str>    NCBI Email account address"
  echo "  --api-key <key>  NCBI API key"
  echo "  --export <str>   Additional environment exports (e.g., 'VAR1=val;VAR2=val')"
  echo "  --dynamic-resources  Dynamically estimate time/mem/CPUs per job"
  echo "  -h               Show this help message"
  exit 1
}

########################################
# Parse command-line flags
########################################
while [[ $# -gt 0 ]]; do
  case $1 in
    -i) ITEMS_FILE="$2"; shift 2 ;;
    -o) COMPLETED_FILE="$2"; CHECKPOINT_FILE="$2"; shift 2 ;;
    -m) MAX_JOBS="$2"; shift 2 ;;
    -p) JOB_NAME_PREFIX="$2"; shift 2 ;;
    -c) CPUS="$2"; shift 2 ;;
    -C) JOB_COMMAND="$2"; shift 2 ;;
    -T) SLURM_TIME="$2"; shift 2 ;;
    -M) SLURM_MEM="$2"; shift 2 ;;
    --api-key) NCBI_API_KEY="$2"; shift 2 ;;
    --email) EMAIL="$2"; shift 2 ;;
    --export) USER_EXPORTS="$2"; shift 2 ;;
    --dynamic-resources) DYNAMIC_RESOURCES=1; shift ;;
    -h) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

declare -A COMPLETED
if [[ -f "$CHECKPOINT_FILE" ]]; then
  while read -r line; do
    line=$(echo "$line" | xargs)
    [[ -n "$line" ]] && COMPLETED["$line"]=1
  done < "$CHECKPOINT_FILE"
fi
########################################
# Apply user exports early
########################################
if [[ -n "$USER_EXPORTS" ]]; then
  IFS=';' read -r -a exports <<< "$USER_EXPORTS"
  for export in "${exports[@]}"; do
    export=$(echo "$export" | xargs)
    [[ -z "$export" ]] && continue
    eval "$export"
  done
fi

# Update dependent variables after user exports
LOG_DIR="${WORKDIR}/logs"
JOBS_DIR="${WORKDIR}/jobs"
TOKEN_FILE="${TOKEN_FILE:-${LOG_DIR}/.job_tokens}"
TOKEN_LOCK_FILE="${TOKEN_LOCK_FILE:-${LOG_DIR}/.job_tokens.lock}"
CHECKPOINT_FILE="${CHECKPOINT_FILE:-${WORKDIR}/completed_accessions.txt}"
CHECKPOINT_LOCK_FILE="${CHECKPOINT_LOCK_FILE:-${LOG_DIR}/checkpoint.lock}"
DEBUG_LOCK_FILE="${DEBUG_LOCK_FILE:-${LOG_DIR}/debug.lock}"
PERIODIC_RESET_LOG="${LOG_DIR}/periodic_reset.log"

########################################
# Setup directories/files and checks
########################################
echo "Setting up directories in $WORKDIR" >&2
for dir in "$JOBS_DIR" "$LOG_DIR" "$(dirname "$CHECKPOINT_FILE")" "$(dirname "$CHECKPOINT_LOCK_FILE")" "$(dirname "$TOKEN_FILE")" "$(dirname "$TOKEN_LOCK_FILE")" "$(dirname "$DEBUG_LOCK_FILE")"; do
  mkdir -p "$dir" || { echo "ERROR: Failed to create $dir" >&2; exit 1; }
done

for file in "$CHECKPOINT_FILE" "$CHECKPOINT_LOCK_FILE" "$TOKEN_FILE" "$DEBUG_LOCK_FILE"; do
  touch "$file" 2>/dev/null || { echo "ERROR: Failed to create $file" >&2; exit 1; }
done

# Debug file locations
{
  echo "[$(date)] Setup:"
  echo "WORKDIR=$WORKDIR"
  echo "LOG_DIR=$LOG_DIR"
  echo "TOKEN_FILE=$TOKEN_FILE"
  echo "TOKEN_LOCK_FILE=$TOKEN_LOCK_FILE"
  echo "CHECKPOINT_FILE=$CHECKPOINT_FILE"
  echo "CHECKPOINT_LOCK_FILE=$CHECKPOINT_LOCK_FILE"
  echo "DEBUG_LOCK_FILE=$DEBUG_LOCK_FILE"
} >> "${LOG_DIR}/setup_debug.log"

# Check for required tools
for cmd in sbatch squeue flock curl python3; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd is required but not found" >&2; exit 1; }
done

# Initialize or adjust token file
if [[ ! -f "$TOKEN_FILE" ]]; then
  # Token file doesn’t exist, initialize it with MAX_JOBS
  (
    flock -x 9
    echo "$MAX_JOBS" > "$TOKEN_FILE"
    echo "[$(date)] Initialized token file with $MAX_JOBS tokens" >> "${LOG_DIR}/token_audit.log"
  ) 9>"$TOKEN_LOCK_FILE"
else
  # Token file exists, check running jobs and adjust if needed
  (
    flock -x 9
    current_tokens=$(cat "$TOKEN_FILE")
    running_jobs=$(squeue --noheader --format '%j' -u "$USER" 2>/dev/null | grep -c "^${JOB_NAME_PREFIX}" || echo "0")
    expected_tokens=$((MAX_JOBS - running_jobs))

    if [[ ! "$current_tokens" =~ ^[0-9]+$ || $current_tokens -ne $expected_tokens ]]; then
      echo "$expected_tokens" > "$TOKEN_FILE"
      echo "[$(date)] Reset token file to $expected_tokens (running jobs: $running_jobs)" >> "${LOG_DIR}/token_audit.log"
    else
      echo "[$(date)] Token file OK: $current_tokens tokens, $running_jobs running jobs" >> "${LOG_DIR}/token_audit.log"
    fi
  ) 9>"$TOKEN_LOCK_FILE"
fi

########################################
# periodic_reset function
########################################
periodic_reset() {
  (
    flock -x 9 || { echo "ERROR: Failed to lock $TOKEN_LOCK_FILE" >&2; exit 1; }
    
    # Read current tokens, default to 0 if invalid
    current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null || echo "0")
    if [[ ! "$current_tokens" =~ ^[0-9]+$ ]]; then
      current_tokens=0
    fi
    
    # Count running and pending jobs
    active_jobs=$(squeue --noheader --format '%j' -u "$USER" -t RUNNING,PENDING 2>/dev/null | grep -c "^${JOB_NAME_PREFIX}" || echo "0")
    
    # Calculate expected tokens
    expected_tokens=$((MAX_JOBS - active_jobs))
    if (( expected_tokens < 0 )); then
      expected_tokens=0  # Prevent negative tokens
    fi
    
    # Update token file if misaligned
    if (( current_tokens != expected_tokens )); then
      echo "$expected_tokens" > "$TOKEN_FILE"
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reset tokens to $expected_tokens (active jobs: $active_jobs)" >> "$PERIODIC_RESET_LOG"
    else
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Tokens OK: $current_tokens (active jobs: $active_jobs)" >> "$PERIODIC_RESET_LOG"
    fi
  ) 9>"$TOKEN_LOCK_FILE" || { echo "ERROR: Failed to reset tokens" >&2; return 1; }
  return 0
}

########################################
# cleanup_successful_jobs function
########################################
cleanup_successful_jobs() {
  (
    flock -n 9 || return 0  # Skip if locked, no big deal
    
    # Process the last 5000 completed items
    tail -n 5000 "$CHECKPOINT_FILE" | while read -r ITEM; do
      ITEM=$(echo "$ITEM" | xargs)
      [[ -z "$ITEM" ]] && continue
      
      # Clean up files, track deletions
      deleted=0
      for file in "${JOBS_DIR}/${ITEM}.sh" "${LOG_DIR}/${ITEM}.out" "${LOG_DIR}/${ITEM}.err" "${LOG_DIR}/${ITEM}.debug.log"; do
        [[ -f "$file" ]] && rm -f "$file" && ((deleted++))
      done
      [[ $deleted -gt 0 ]] && echo "[$(date)] Cleaned $deleted files for $ITEM" >> "${LOG_DIR}/cleanup.log"
    done
    
    exit 0
  ) 9>"$CHECKPOINT_LOCK_FILE" || return 0
  return 0
}

########################################
# estimate_resources function
########################################
estimate_resources() {
  local ITEM="$1"
  local size_bytes=0
  local time="10:00:00"  # Default: 10 hours
  local mem="8G"         # Default: 8 GB
  local cpus="4"         # Default: 4 CPUs
  local prefix="${ITEM:0:3}"

  # Determine provider
  case "$prefix" in
    SRR|SRX|SRS|SRP)
      # SRA: Single attempt, fallback to default if fails
      size_bytes=$(esummary -db sra -id "$ITEM" 2>/dev/null | xtract -pattern DocumentSummary -element Size 2>/dev/null || echo "0")
      ;;
    ERR|ERX|ERS|ERP|DRR|DRX|DRS|DRP)
      # ENA: Single attempt, sum fastq_bytes
      size_bytes=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ITEM}&result=read_run&fields=fastq_bytes" \
        | awk 'NR>1 {split($1, a, ","); for (i in a) s+=a[i]} END {print s+0}' 2>/dev/null || echo "0")
      ;;
    *)
      echo "[$(date)] Unknown provider for $ITEM" >> "${LOG_DIR}/resource_errors.log"
      echo "$time $mem $cpus"
      return
      ;;
  esac

  # Validate size, log if fetch failed
  if [[ ! "$size_bytes" =~ ^[0-9]+$ || $size_bytes -eq 0 ]]; then
    echo "[$(date)] Failed to fetch size for $ITEM, using defaults" >> "${LOG_DIR}/resource_errors.log"
    size_bytes=0
  fi

  # Assign resources based on size
  if (( size_bytes > 20000000000 )); then      # >20 GB
    time="12:00:00"; mem="12G"; cpus="8"
  elif (( size_bytes > 5000000000 )); then    # >5 GB
    time="04:00:00"; mem="8G"; cpus="6"
  elif (( size_bytes > 1000000000 )); then    # >1 GB
    time="02:00:00"; mem="6G"; cpus="4"
  else                                        # ≤1 GB or unknown
    time="01:00:00"; mem="4G"; cpus="2"
  fi

  echo "$time $mem $cpus"
}

########################################
# submit_job function
########################################
submit_job() {
  local ITEM="$1"
  local JOB_SCRIPT="${JOBS_DIR}/${ITEM}.sh"

  ITEM=$(echo "$ITEM" | xargs)
  [[ -z "$ITEM" ]] && return

  ENV_EXPORTS="export ITEM=\"${ITEM}\""
  ENV_EXPORTS+=$'\n'"export ACCESSION=\"${ITEM}\""
  ENV_EXPORTS+=$'\n'"export WORKDIR=\"${WORKDIR}\""
  ENV_EXPORTS+=$'\n'"export LOG_DIR=\"${LOG_DIR}\""
  ENV_EXPORTS+=$'\n'"export CHECKPOINT_FILE=\"${CHECKPOINT_FILE}\""
  ENV_EXPORTS+=$'\n'"export CHECKPOINT_LOCK_FILE=\"${CHECKPOINT_LOCK_FILE}\""
  ENV_EXPORTS+=$'\n'"export TOKEN_FILE=\"${TOKEN_FILE}\""
  ENV_EXPORTS+=$'\n'"export TOKEN_LOCK_FILE=\"${TOKEN_LOCK_FILE}\""
  ENV_EXPORTS+=$'\n'"export DEBUG_LOCK=\"${DEBUG_LOCK_FILE}\""
  ENV_EXPORTS+=$'\n'"export COMBINED_METADATA=\"${COMBINED_METADATA:-${WORKDIR}/combined_metadata.tsv}\""
  [[ -n "$NCBI_API_KEY" ]] && ENV_EXPORTS+=$'\n'"export NCBI_API_KEY=\"${NCBI_API_KEY}\""
  [[ -n "$EMAIL" ]] && ENV_EXPORTS+=$'\n'"export EMAIL=\"${EMAIL}\""
  [[ -n "$USER_EXPORTS" ]] && for export in "${exports[@]}"; do ENV_EXPORTS+=$'\n'"export $export"; done

  cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=${JOB_NAME_PREFIX}${ITEM}
#SBATCH --output=${LOG_DIR}/${ITEM}.out
#SBATCH --error=${LOG_DIR}/${ITEM}.err
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --ntasks=1
#SBATCH --account=one

set -euo pipefail

# Environment exports
${ENV_EXPORTS}

# Debug environment variables
echo "Running job for \$ITEM" >> "${LOG_DIR}/${ITEM}.debug.log"
for var in WORKDIR CHECKPOINT_FILE CHECKPOINT_LOCK_FILE TOKEN_FILE TOKEN_LOCK_FILE DEBUG_LOCK; do
  printf '%s=%s\n' "\$var" "\${!var}" >> "${LOG_DIR}/${ITEM}.debug.log"
done

# Function to release token using flock timeout
release_token() {
  # Open file descriptor 9 for the lock file
  exec 9>"$TOKEN_LOCK_FILE" || {
    echo "ERROR: Unable to open lock file $TOKEN_LOCK_FILE" >> "${LOG_DIR}/${ITEM}.err"
    return 1
  }

  # Attempt to acquire an exclusive lock, waiting up to 10 seconds
  if ! flock -x -w 10 9; then
    echo "ERROR: Failed to acquire lock for $ITEM" >> "${LOG_DIR}/${ITEM}.err"
    return 1
  fi

  # Read the current token value, defaulting to 0 if missing
  current=$(<"$TOKEN_FILE" 2>/dev/null || echo "0")

  # Validate that the token is a non-negative integer
  if ! [[ "$current" =~ ^[0-9]+$ ]]; then
    echo "Invalid token count: $current" >> "${LOG_DIR}/token_errors.log"
    flock -u 9
    return 1
  fi

  # Increment the token and write the new value back to the file
  new=$((current + 1))
  echo "$new" > "$TOKEN_FILE"

  # Log the token release event with a timestamp
  echo "[$(date)] RELEASED TOKEN FOR $ITEM (NOW $new)" >> "${LOG_DIR}/token_audit.log"

  # Release the lock
  flock -u 9
  return 0
}

# Ensure that the token is released when the script exits
trap 'release_token' EXIT

# Run the command with retry
EXIT_CODE=0
set +e
${JOB_COMMAND}
EXIT_CODE=\$?
set -e
echo "Final exit code: \$EXIT_CODE" >> "${LOG_DIR}/${ITEM}.debug.log"
exit \$EXIT_CODE
EOF

  chmod +x "$JOB_SCRIPT" || { echo "ERROR: Failed to make $JOB_SCRIPT executable" >&2; return 1; }
  sbatch "$JOB_SCRIPT" && { echo "Submitted job for item: $ITEM" >&2; return 0; }
  return 1
}

########################################
# Main logic
########################################
CLEANUP_COUNTER=0

# Trap for cleanup on exit
trap 'cleanup_successful_jobs; echo "Script terminated" >&2' EXIT

while IFS= read -r ITEM || [[ -n "$ITEM" ]]; do
  CURRENT_TIME=$(date +%s)
  ELAPSED=$(( CURRENT_TIME - START_TIME ))
  if (( ELAPSED > MAX_RUNTIME )); then
    echo "Time nearly up (47:59). Exiting for resubmission..." >&2
    exit 42
  fi

  ITEM=$(echo "$ITEM" | xargs)
  [[ -z "$ITEM" ]] && continue
  [[ -n "${COMPLETED[$ITEM]}" ]] && continue

  RETRY_COUNT=0
  MAX_RETRIES=3600
    while (( RETRY_COUNT < MAX_RETRIES )); do
    (
      flock -x 9 || { echo "Lock starvation detected" >> "${LOG_DIR}/lock_errors.log"; exit 99; }
      current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null || echo "0")
      if ! [[ "$current_tokens" =~ ^[0-9]+$ ]]; then
        echo "Invalid token count: $current_tokens" >> "${LOG_DIR}/token_errors.log"
        echo "$MAX_JOBS" > "$TOKEN_FILE"
        current_tokens=$MAX_JOBS
      fi
      if (( current_tokens > 0 )); then
        echo $((current_tokens - 1)) > "$TOKEN_FILE"
        echo "[$(date)] ACQUIRED TOKEN FOR $ITEM (now $((current_tokens - 1)))" >> "${LOG_DIR}/token_audit.log"
        exit 0
      fi
      exit 1
    ) 9>"$TOKEN_LOCK_FILE"

    case $? in
      0) break ;;
      1)
        sleep 1
        CURRENT_TIME=$(date +%s)
        if (( CURRENT_TIME - LAST_RESET_TIME >= RESET_INTERVAL )); then
          if periodic_reset; then
            LAST_RESET_TIME=$CURRENT_TIME
          else
            echo "Periodic reset failed" >&2
            exit 1
          fi
        fi
        ((RETRY_COUNT++))
        ;;
      99)
        sleep 2
        ((RETRY_COUNT++))
        ;;
    esac
  done
  if (( RETRY_COUNT >= MAX_RETRIES )); then
    echo "ERROR: Failed to acquire token for $ITEM after $MAX_RETRIES attempts" >&2
    exit 1
  fi

  [[ "$DYNAMIC_RESOURCES" == "1" ]] && read SLURM_TIME SLURM_MEM CPUS < <(estimate_resources "$ITEM")
  submit_job "$ITEM" || { echo "ERROR: Job submission failed for $ITEM" >&2; exit 1; }

  if (( ++CLEANUP_COUNTER >= 100 )); then
    cleanup_successful_jobs &
    CLEANUP_COUNTER=0
  fi
done < "$ITEMS_FILE"

cleanup_successful_jobs

echo "🎉 All jobs completed!" >&2
