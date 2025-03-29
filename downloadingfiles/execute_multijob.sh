#!/bin/bash
#
# A general-purpose Slurm batch “driver” script that:
#   1) Reads tasks (one per line) from an input file
#   2) Submits each task to Slurm with concurrency-limiting tokens
#   3) Records completed tasks in a checkpoint file
#   4) Periodically cleans up leftover logs/job scripts
#   5) Exports environment variables (like NCBI_API_KEY) to each job
#

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
MAX_JOBS=20
JOB_NAME_PREFIX="concurrent_"
JOB_COMMAND="echo 'Processing \$ITEM'"
SLURM_TIME="02:00:00"
SLURM_MEM="8G"
CPUS="4"
THIS_JOB_ID="${SLURM_JOB_ID}"
ORIGINAL_COMMAND="sbatch $0 $@"
START_TIME="$(date +%s)"
MAX_RUNTIME=$(( 47 * 3600 - 60 ))  # 47 hours minus 1 minute

NCBI_API_KEY=""
EMAIL=""
USER_EXPORTS=""
DYNAMIC_RESOURCES=0  

WORKDIR="$(pwd)"
LOG_DIR="${WORKDIR}/logs"
JOBS_DIR="${WORKDIR}/jobs"
TOKEN_FILE="${LOG_DIR}/.job_tokens"
TOKEN_LOCK_FILE="${LOG_DIR}/.job_tokens.lock"
CHECKPOINT_LOCK_FILE="${LOG_DIR}/checkpoint.lock"
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
  echo "  --email <str>    NCBI Email account address (your-email@email.edu)"
  echo "  --api-key <key>  NCBI API key or any other single-value key"
  echo "  --export <str>   Additional environment exports (e.g. 'VAR1=val VAR2=val')"
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
    -o) COMPLETED_FILE="$2"; shift 2 ;;
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

########################################
# Setup directories/files and dependency checks
########################################
mkdir -p "$JOBS_DIR" "$LOG_DIR"
touch "$COMPLETED_FILE" "$CHECKPOINT_LOCK_FILE"

# Check for required tools
for cmd in sbatch squeue flock curl; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: $cmd is required but not found"; exit 1; }
done

# Initialize token file if missing or invalid
if ! [[ -f "$TOKEN_FILE" && "$(cat "$TOKEN_FILE")" =~ ^[0-9]+$ ]]; then
  echo "$MAX_JOBS" > "$TOKEN_FILE"
fi

########################################
# periodic_reset function
########################################
periodic_reset() {
  {
    flock -x 9 || return 1
    current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null)
    if ! [[ "$current_tokens" =~ ^[0-9]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $TOKEN_FILE invalid: '$current_tokens'. Resetting to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
      echo "$MAX_JOBS" > "$TOKEN_FILE"
      return 0
    fi
  } 9<>"$TOKEN_LOCK_FILE"

  RUNNING_JOBS=$(squeue --noheader --format '%j' 2>/dev/null | grep -c "^${JOB_NAME_PREFIX}" || echo "0")
  tokens_out=$((MAX_JOBS - current_tokens))

  if (( tokens_out < 0 )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: tokens_out < 0 => reset to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
    {
      flock -x 9 || return 1
      echo "$MAX_JOBS" > "$TOKEN_FILE"
    } 9<>"$TOKEN_LOCK_FILE"
    return 0
  fi

  if (( tokens_out > RUNNING_JOBS )); then
    mismatch=$(( tokens_out - RUNNING_JOBS ))
    new_tokens=$(( current_tokens + mismatch ))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] MISMATCH: tokens_out=$tokens_out, running=$RUNNING_JOBS. Restoring $mismatch tokens." | tee -a "$PERIODIC_RESET_LOG"
    {
      flock -x 9 || return 1
      echo "$new_tokens" > "$TOKEN_FILE"
    } 9<>"$TOKEN_LOCK_FILE"
  else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK: tokens_out=$tokens_out, running=$RUNNING_JOBS, current_tokens=$current_tokens" >> "$PERIODIC_RESET_LOG"
  fi
}

########################################
# cleanup_successful_jobs function
########################################
cleanup_successful_jobs() {
  flock -n "$CHECKPOINT_LOCK_FILE" || { echo "Cleanup skipped: lock held"; return; }
  tail -n 100 "$COMPLETED_FILE" | while read -r ITEM; do
    ITEM=$(echo "$ITEM" | xargs)
    [[ -z "$ITEM" ]] && continue

    SCRIPT_PATH="${JOBS_DIR}/${ITEM}.sh"
    [[ -f "$SCRIPT_PATH" ]] && rm -f "$SCRIPT_PATH" && echo "Removed job script for $ITEM"

    LOG_OUT="${LOG_DIR}/${ITEM}.out"
    LOG_ERR="${LOG_DIR}/${ITEM}.err"
    [[ -f "$LOG_OUT" ]] && rm -f "$LOG_OUT" && echo "Removed log file $LOG_OUT"
    [[ -f "$LOG_ERR" ]] && rm -f "$LOG_ERR" && echo "Removed log file $LOG_ERR"
  done
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

  if [[ "$provider" == "sra" ]]; then
    size_bytes=$(esummary -db sra -id "$ITEM" 2>/dev/null | xtract -pattern DocumentSummary -element Size 2>/dev/null)
  else
    size_bytes=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ITEM}&result=read_run&fields=fastq_bytes" \
      | tail -n +2 | tr ',' '\n' | awk '{s+=$1} END {print s}' 2>/dev/null)
  fi

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
  ENV_EXPORTS+=$'\n'"export LOG_DIR=\"${LOG_DIR}\""
  ENV_EXPORTS+=$'\n'"export CHECKPOINT_LOCK_FILE=\"${CHECKPOINT_LOCK_FILE}\""
  ENV_EXPORTS+=$'\n'"export TOKEN_LOCK_FILE=\"${TOKEN_LOCK_FILE}\""
  ENV_EXPORTS+=$'\n'"export COMPLETED_FILE=\"${COMPLETED_FILE}\""

  [[ -n "$NCBI_API_KEY" ]] && ENV_EXPORTS+=$'\n'"export NCBI_API_KEY=\"${NCBI_API_KEY}\""
  [[ -n "$EMAIL" ]] && ENV_EXPORTS+=$'\n'"export EMAIL=\"${EMAIL}\""
  [[ -n "$USER_EXPORTS" ]] && for kv in $USER_EXPORTS; do ENV_EXPORTS+=$'\n'"export $kv"; done

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

# Run the command
EXIT_CODE=0
set +e
${JOB_COMMAND}
EXIT_CODE=$?
set -e
echo "Exit code: $EXIT_CODE" >> "${LOG_DIR}/${ITEM}.debug.log"

# Release token
{
  flock -x 9
  current_tokens=\$(cat "\$TOKEN_FILE")
  echo \$((current_tokens + 1)) > "\$TOKEN_FILE"
  echo "[$(date)] RELEASED TOKEN FOR \$ITEM (now \$((current_tokens + 1)))" >> "\${LOG_DIR}/token_audit.log"
} 9>"\$TOKEN_LOCK_FILE"

# Checkpoint only on success
if [[ \$EXIT_CODE -eq 0 ]]; then
  {
    flock 200
    echo "\${ITEM}" >> "\${COMPLETED_FILE}"
  } 200>"\${CHECKPOINT_LOCK_FILE}"
else
  echo "Job failed, not checkpointing" >> "\${LOG_DIR}/${ITEM}.debug.log"
  exit \$EXIT_CODE
fi
EOF

  chmod +x "$JOB_SCRIPT"
  sbatch "$JOB_SCRIPT"
  echo "Submitted job for item: $ITEM"
}

########################################
# Main logic
########################################
CLEANUP_COUNTER=0

while IFS= read -r ITEM; do
  CURRENT_TIME=$(date +%s)
  ELAPSED=$(( CURRENT_TIME - START_TIME ))
  if (( ELAPSED > MAX_RUNTIME )); then
    echo "Time nearly up (47:59). Exiting for resubmission..." >&2
    exit 42  # Signal wrapper script to resubmit
  fi

  grep -Fxq "$ITEM" "$COMPLETED_FILE" && continue
  ITEM=$(echo "$ITEM" | xargs)
  [[ -z "$ITEM" ]] && continue

  # Acquire token with timeout
  RETRY_COUNT=0
  MAX_RETRIES=3600
  while (( RETRY_COUNT < MAX_RETRIES )); do
    (
      flock -x 9 || exit 99
      current_tokens=$(< "$TOKEN_FILE")
      if (( current_tokens > 0 )); then
        echo $((current_tokens - 1)) > "$TOKEN_FILE"
        echo "[$(date)] ACQUIRED TOKEN FOR $ITEM (now $((current_tokens - 1)))" >> "${LOG_DIR}/token_audit.log"
        exit 0
      fi
      exit 1
    ) 9>"$TOKEN_LOCK_FILE"

    case $? in
      0) break ;;
      1) sleep 1; periodic_reset; ((RETRY_COUNT++)) ;;
      99) echo "Lock starvation detected" >> "${LOG_DIR}/lock_errors.log"; sleep 2; ((RETRY_COUNT++)) ;;
    esac
  done
  if (( RETRY_COUNT >= MAX_RETRIES )); then
    echo "ERROR: Failed to acquire token for $ITEM after $MAX_RETRIES attempts" >&2
    exit 1
  fi

  [[ "$DYNAMIC_RESOURCES" == "1" ]] && read SLURM_TIME SLURM_MEM CPUS < <(estimate_resources "$ITEM")
  submit_job "$ITEM"

  if (( ++CLEANUP_COUNTER >= 20 )); then
    cleanup_successful_jobs &
    CLEANUP_COUNTER=0
  fi
done < "$ITEMS_FILE"

# Final cleanup
cleanup_successful_jobs

# Wait for all jobs to finish
while [[ $(flock -x "$TOKEN_LOCK_FILE" -c "cat $TOKEN_FILE") -lt $MAX_JOBS ]]; do
  periodic_reset
  sleep 5
done

echo "🎉 All jobs completed!"
