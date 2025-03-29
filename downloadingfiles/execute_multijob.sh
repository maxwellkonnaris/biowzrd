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

# Initialize token file with retry
for attempt in {1..3}; do
  if ! [[ -f "$TOKEN_FILE" && "$(cat "$TOKEN_FILE")" =~ ^[0-9]+$ ]]; then
    (
      flock -x 9 || { echo "ERROR: Failed to lock $TOKEN_LOCK_FILE" >&2; exit 1; }
      echo "$MAX_JOBS" > "$TOKEN_FILE" || { echo "ERROR: Failed to initialize $TOKEN_FILE" >&2; exit 1; }
    ) 9>"$TOKEN_LOCK_FILE" && break
  fi
  [[ $attempt -lt 3 ]] && sleep 2
done || { echo "ERROR: Failed to initialize tokens after retries" >&2; exit 1; }

########################################
# periodic_reset function
########################################
periodic_reset() {
  local attempt
  for attempt in {1..3}; do
    (
      flock -x 9 || { echo "ERROR: Failed to lock $TOKEN_LOCK_FILE (attempt $attempt)" >> "${LOG_DIR}/lock_errors.log"; exit 1; }
      current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null || echo "0")
      if ! [[ "$current_tokens" =~ ^[0-9]+$ ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $TOKEN_FILE invalid: '$current_tokens'. Resetting to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
        echo "$MAX_JOBS" > "$TOKEN_FILE"
        exit 0
      fi

      RUNNING_JOBS=$(squeue --noheader --format '%j' -u "$USER" 2>/dev/null | grep -c "^${JOB_NAME_PREFIX}" || echo "0")
      tokens_out=$((MAX_JOBS - current_tokens))

      if (( tokens_out < 0 )); then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: tokens_out < 0 => reset to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
        echo "$MAX_JOBS" > "$TOKEN_FILE"
        exit 0
      fi

      if (( tokens_out > RUNNING_JOBS )); then
        mismatch=$(( tokens_out - RUNNING_JOBS ))
        new_tokens=$(( current_tokens + mismatch ))
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] MISMATCH: tokens_out=$tokens_out, running=$RUNNING_JOBS. Restoring $mismatch tokens." | tee -a "$PERIODIC_RESET_LOG"
        echo "$new_tokens" > "$TOKEN_FILE"
      elif (( RUNNING_JOBS > MAX_JOBS )); then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: running=$RUNNING_JOBS exceeds MAX_JOBS=$MAX_JOBS. Limiting to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
        echo "0" > "$TOKEN_FILE"  # Force no new jobs until running drops
      else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK: tokens_out=$tokens_out, running=$RUNNING_JOBS, current_tokens=$current_tokens" >> "$PERIODIC_RESET_LOG"
      fi
      exit 0
    ) 9>"$TOKEN_LOCK_FILE" && return 0
    sleep $((attempt * 2))
  done
  echo "ERROR: Failed to reset tokens after retries" >&2
  return 1
}

########################################
# cleanup_successful_jobs function
########################################
cleanup_successful_jobs() {
  local attempt
  for attempt in {1..3}; do
    (
      flock -n 9 || { echo "Cleanup skipped: lock held (attempt $attempt)" >> "${LOG_DIR}/lock_errors.log"; exit 1; }
      tail -n 100 "$CHECKPOINT_FILE" | while read -r ITEM; do
        ITEM=$(echo "$ITEM" | xargs)
        [[ -z "$ITEM" ]] && continue

        SCRIPT_PATH="${JOBS_DIR}/${ITEM}.sh"
        [[ -f "$SCRIPT_PATH" ]] && rm -f "$SCRIPT_PATH" && echo "Removed job script for $ITEM" >> "${LOG_DIR}/cleanup.log"

        LOG_OUT="${LOG_DIR}/${ITEM}.out"
        LOG_ERR="${LOG_DIR}/${ITEM}.err"
        [[ -f "$LOG_OUT" ]] && rm -f "$LOG_OUT" && echo "Removed log file $LOG_OUT" >> "${LOG_DIR}/cleanup.log"
        [[ -f "$LOG_ERR" ]] && rm -f "$LOG_ERR" && echo "Removed log file $LOG_ERR" >> "${LOG_DIR}/cleanup.log"
      done
      exit 0
    ) 9>"$CHECKPOINT_LOCK_FILE" && return 0
    sleep $((attempt * 2))
  done
  echo "WARNING: Cleanup failed after retries" >&2
  return 1
}

########################################
# estimate_resources function
########################################
estimate_resources() {
  local ITEM="$1"
  local provider=""
  local size_bytes=0
  local time="01:00:00"
  local mem="1G"
  local cpus="1"

  local prefix="${ITEM:0:3}"
  case "$prefix" in
    SRR|SRX|SRS|SRP) provider="sra" ;;
    ERR|ERX|ERS|ERP|DRR|DRX|DRS|DRP) provider="ena" ;;
    *) echo "Unknown provider for accession $ITEM" >> "${LOG_DIR}/resource_errors.log"; echo "$time $mem $cpus"; return ;;
  esac

  for attempt in {1..3}; do
    if [[ "$provider" == "sra" ]]; then
      size_bytes=$(esummary -db sra -id "$ITEM" 2>/dev/null | xtract -pattern DocumentSummary -element Size 2>/dev/null)
    else
      size_bytes=$(curl -s --retry 3 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ITEM}&result=read_run&fields=fastq_bytes" \
        | tail -n +2 | tr ',' '\n' | awk '{s+=$1} END {print s}' 2>/dev/null)
    fi
    [[ -n "$size_bytes" && "$size_bytes" =~ ^[0-9]+$ ]] && break
    echo "Retry $attempt: Failed to estimate size for $ITEM" >> "${LOG_DIR}/resource_errors.log"
    sleep $((attempt * 2))
  done

  [[ -z "$size_bytes" || ! "$size_bytes" =~ ^[0-9]+$ ]] && { echo "Failed to estimate size for $ITEM" >> "${LOG_DIR}/resource_errors.log"; size_bytes=0; }

  if (( size_bytes > 20000000000 )); then      # >20 GB
    time="12:00:00"; mem="10G"; cpus="10"
  elif (( size_bytes > 5000000000 )); then    # >5 GB
    time="03:00:00"; mem="6G"; cpus="6"
  else                                        # Small
    time="01:00:00"; mem="4G"; cpus="4"
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

# Function to release token
release_token() {
  for attempt in {1..5}; do
    (
      flock -x 9 || exit 1
      current=\$(cat "\$TOKEN_FILE" 2>/dev/null || echo "0")
      if [[ "\$current" =~ ^[0-9]+$ ]]; then
        new=\$((current + 1))
        echo "\$new" > "\$TOKEN_FILE"
        echo "[$(date)] RELEASED TOKEN FOR \$ITEM (NOW \$new)" >> "${LOG_DIR}/token_audit.log"
        exit 0
      fi
      echo "Invalid token count: \$current" >> "${LOG_DIR}/token_errors.log"
      exit 1
    ) 9>"\$TOKEN_LOCK_FILE" && return 0
    sleep \$((attempt * 2))
  done
  echo "ERROR: Failed to release token for \$ITEM after 5 attempts" >> "${LOG_DIR}/${ITEM}.err"
}

# Trap to ensure token release on exit
trap 'release_token' EXIT

# Run the command with retry
EXIT_CODE=0
for attempt in {1..3}; do
  set +e
  ${JOB_COMMAND}
  EXIT_CODE=\$?
  set -e
  [[ \$EXIT_CODE -eq 0 ]] && break
  echo "Attempt \$attempt failed with exit code \$EXIT_CODE. Retrying..." >> "${LOG_DIR}/${ITEM}.debug.log"
  sleep \$((attempt * 5))
done
echo "Final exit code: \$EXIT_CODE" >> "${LOG_DIR}/${ITEM}.debug.log"
exit \$EXIT_CODE
EOF

  chmod +x "$JOB_SCRIPT" || { echo "ERROR: Failed to make $JOB_SCRIPT executable" >&2; return 1; }
  for attempt in {1..3}; do
    sbatch "$JOB_SCRIPT" && { echo "Submitted job for item: $ITEM" >&2; return 0; }
    echo "ERROR: Failed to submit $JOB_SCRIPT (attempt $attempt)" >&2
    sleep $((attempt * 2))
  done
  echo "ERROR: Failed to submit $JOB_SCRIPT after retries" >&2
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

  if (( ++CLEANUP_COUNTER >= 20 )); then
    cleanup_successful_jobs &
    CLEANUP_COUNTER=0
  fi
done < "$ITEMS_FILE"

cleanup_successful_jobs

# Wait for all tokens to return with timeout
TIMEOUT=$((60 * 5))  # 5 minutes
START_WAIT=$(date +%s)
while true; do
  CURRENT_WAIT=$(date +%s)
  ELAPSED_WAIT=$((CURRENT_WAIT - START_WAIT))
  if (( ELAPSED_WAIT > TIMEOUT )); then
    echo "ERROR: Timeout waiting for tokens to return (current_tokens=$(cat "$TOKEN_FILE"))" >&2
    exit 1
  fi
  (
    flock -x 9 || exit 1
    current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null || echo "0")
    [[ "$current_tokens" =~ ^[0-9]+$ ]] || { echo "$MAX_JOBS" > "$TOKEN_FILE"; current_tokens=$MAX_JOBS; }
    [[ "$current_tokens" -eq "$MAX_JOBS" ]] && exit 0
    exit 1
  ) 9>"$TOKEN_LOCK_FILE" && break
  periodic_reset || { echo "Periodic reset failed in final wait" >&2; exit 1; }
  sleep 5
done

echo "🎉 All jobs completed!" >&2
