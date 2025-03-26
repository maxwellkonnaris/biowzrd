#!/bin/bash
#
# A general-purpose Slurm batch “driver” script that:
#   1) Reads tasks (one per line) from an input file
#   2) Submits each task to Slurm with concurrency-limiting tokens
#   3) Records completed tasks in a checkpoint file
#   4) Periodically cleans up leftover logs/job scripts
#   5) Exports environment variables (like NCBI_API_KEY) to each job
#
# Example usage:
#   sbatch execute_slurmmultijob.sh \
#       -i run_accessions.txt \
#       -o completed_accessions.txt \
#       -m 20 \
#       -p fastq_ \
#       -C "python /path/to/download_accession.py" \
#       -T "04:00:00" \
#       -M "8G" \
#       --email "XXXXX" \
#       --api-key "XXXXXX" \
#       --export "WORKDIR=$(pwd)"


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
MAX_RUNTIME=$(( 47 * 3600 -60 ))


# Let user optionally define an API key, plus arbitrary "export" lines
NCBI_API_KEY=""
EMAIL=""
USER_EXPORTS=""
DYNAMIC_RESOURCES=0  

WORKDIR="$(pwd)"
LOG_DIR="${WORKDIR}/logs"
JOBS_DIR="${WORKDIR}/jobs"
#MET_DIR="${WORKDIR}/metadata"

TOKEN_FILE="${WORKDIR}/.job_tokens"
TOKEN_LOCK_FILE="${WORKDIR}/.job_tokens.lock"
CHECKPOINT_LOCK_FILE="${WORKDIR}/checkpoint.lock"
PERIODIC_RESET_LOG="${LOG_DIR}/periodic_reset.log"
CLEANUP_COUNTER=0

########################################
# usage function
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
    -i)
      ITEMS_FILE="$2"
      shift 2
      ;;
    -o)
      COMPLETED_FILE="$2"
      shift 2
      ;;
    -m)
      MAX_JOBS="$2"
      shift 2
      ;;
    -p)
      JOB_NAME_PREFIX="$2"
      shift 2
      ;;
    -c)
      CPUS="$2"
      shift 2
      ;;
    -C)
      JOB_COMMAND="$2"
      shift 2
      ;;
    -T)
      SLURM_TIME="$2"
      shift 2
      ;;
    -M)
      SLURM_MEM="$2"
      shift 2
      ;;
    --api-key)
      NCBI_API_KEY="$2"
      shift 2
      ;;
    --email)
      EMAIL="$2"
      shift 2
      ;;
    --export)
      USER_EXPORTS="$2"
      shift 2
      ;;
    --dynamic-resources)
      DYNAMIC_RESOURCES=1
      shift
      ;;
    -h)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

########################################
# Setup directories/files
########################################
mkdir -p "$JOBS_DIR" "$LOG_DIR" #"$MET_DIR"
touch "$COMPLETED_FILE"
touch "$CHECKPOINT_LOCK_FILE" 

# Initialize the token file if missing
if ! [[ -f "$TOKEN_FILE" && "$(cat "$TOKEN_FILE")" =~ ^[0-9]+$ ]]; then
  echo "$MAX_JOBS" > "$TOKEN_FILE"
fi


########################################
# periodic_reset function
########################################
periodic_reset() {
  {
    flock -x 9 || return 1
    if [[ -s "$TOKEN_FILE" ]]; then
      current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null)
    else
      current_tokens="0"
    fi
    if ! [[ "$current_tokens" =~ ^[0-9]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $TOKEN_FILE content invalid: '$current_tokens'. Resetting to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
      echo "$MAX_JOBS" > "$TOKEN_FILE"
      return 0
    fi
  } 9<>"$TOKEN_LOCK_FILE"

  local RUNNING_JOBS
  RUNNING_JOBS=$(squeue --noheader --format '%j' | grep -c "^${JOB_NAME_PREFIX}")
  local tokens_out=$((MAX_JOBS - current_tokens))

  if (( tokens_out < 0 )); then
    # corrupted token file => forcibly reset
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: tokens_out < 0 => reset to $MAX_JOBS." | tee -a "$PERIODIC_RESET_LOG"
    {
      flock -x 9 || return 1
      echo "$MAX_JOBS" > "$TOKEN_FILE"
    } 9<>"$TOKEN_LOCK_FILE"
    return 0
  fi

  if (( tokens_out > RUNNING_JOBS )); then
    # Mismatch => restore tokens
    local mismatch=$(( tokens_out - RUNNING_JOBS ))
    local new_tokens=$(( current_tokens + mismatch ))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] MISMATCH DETECTED: tokens_out=$tokens_out, running=$RUNNING_JOBS. Restoring $mismatch tokens." | tee -a "$PERIODIC_RESET_LOG"
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
  flock -x "$CHECKPOINT_LOCK_FILE"
  tail -n 100 "$COMPLETED_FILE" | while read -r ITEM; do
    ITEM=$(echo "$ITEM" | xargs)
    [[ -z "$ITEM" ]] && continue

    local SCRIPT_PATH="${JOBS_DIR}/${ITEM}.sh"
    if [[ -f "$SCRIPT_PATH" ]]; then
      rm -f "$SCRIPT_PATH"
      echo "Removed job script for $ITEM"
    fi

    local LOG_OUT="${LOG_DIR}/${ITEM}.out"
    local LOG_ERR="${LOG_DIR}/${ITEM}.err"
    [[ -f "$LOG_OUT" ]] && rm -f "$LOG_OUT" && echo "Removed log file $LOG_OUT"
    [[ -f "$LOG_ERR" ]] && rm -f "$LOG_ERR" && echo "Removed log file $LOG_ERR"
  done
  flock -u "$CHECKPOINT_LOCK_FILE"
}

estimate_resources() {
  local ITEM="$1"
  local provider=""
  local size_bytes=0
  local time="01:00:00"
  local mem="1G"
  local cpus="1"

  local prefix="${ITEM:0:3}"
  case "$prefix" in
    SRR|SRX|SRS|SRP)
      provider="sra"
      ;;
    ERR|ERX|ERS|ERP|DRR|DRX|DRS|DRP)
      provider="ena"
      ;;
    *)
      echo "Unknown provider for accession $ITEM" >&2
      echo "$time $mem $cpus"
      return
      ;;
  esac

  # Estimate size
  if [[ "$provider" == "sra" ]]; then
    size_bytes=$(esummary -db sra -id "$ITEM" | xtract -pattern DocumentSummary -element Size 2>/dev/null)
  else
    size_bytes=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ITEM}&result=read_run&fields=fastq_bytes" \
      | tail -n +2 | tr ',' '\n' | awk '{s+=$1} END {print s}')
  fi

  # Default to 0 if unset or failed
  [[ -z "$size_bytes" || ! "$size_bytes" =~ ^[0-9]+$ ]] && size_bytes=0

  # Set thresholds (tune as needed)
if (( size_bytes > 20000000000 )); then      # >20 GB
    time="12:00:00"; mem="8G"; cpus="8"
  elif (( size_bytes > 5000000000 )); then       # >5 GB
    time="03:00:00"; mem="6G"; cpus="4"
  else                                           # Small
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

  # Build environment exports with proper newlines
  ENV_EXPORTS="export ITEM=\"${ITEM}\""
  ENV_EXPORTS+=$'\n'"export ACCESSION=\"${ITEM}\""


  if [[ -n "$NCBI_API_KEY" ]]; then
    ENV_EXPORTS+=$'\n'"export NCBI_API_KEY=\"${NCBI_API_KEY}\""
  fi

  if [[ -n "$EMAIL" ]]; then
    ENV_EXPORTS+=$'\n'"export EMAIL=\"${EMAIL}\""
  fi

  if [[ -n "$USER_EXPORTS" ]]; then
    for kv in $USER_EXPORTS; do
      ENV_EXPORTS+=$'\n'"export $kv"
    done
  fi

  # Create the job script
  cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=${JOB_NAME_PREFIX}${ITEM}
#SBATCH --output=${LOG_DIR}/${ITEM}.out
#SBATCH --error=${LOG_DIR}/${ITEM}.err
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --ntasks=1

set -euo pipefail

# Environment exports
${ENV_EXPORTS}

# The user-specified command:
${JOB_COMMAND}

# Run the command and capture exit status
EXIT_CODE=0
set +e
${JOB_COMMAND}
EXIT_CODE=$?
set -e
echo "Exit code: \$STATUS" >> "${LOG_DIR}/${ITEM}.debug.log"

# Only write to checkpoint if the command succeeded
if [[ \$STATUS -eq 0 ]]; then
  {
    flock 200
    echo "\${ITEM}" >> "${COMPLETED_FILE}"
  } 200>"${CHECKPOINT_LOCK_FILE}"
else
  echo "Job failed, not checkpointing" >> "${LOG_DIR}/${ITEM}.debug.log"
  exit \$STATUS
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

# Read each line of the items file
while IFS= read -r ITEM; do

  # Check run time
  CURRENT_TIME=$(date +%s)
  ELAPSED=$(( CURRENT_TIME - START_TIME ))
  if (( ELAPSED > MAX_RUNTIME )); then
    echo "Time nearly up (47:59). Resubmitting job..."
    echo "$ORIGINAL_COMMAND"
    $ORIGINAL_COMMAND
    exit 0
  fi
  
  # Skip if already completed
  grep -Fxq "$ITEM" "$COMPLETED_FILE" && continue
  ITEM=$(echo "$ITEM" | xargs)
  [[ -z "$ITEM" ]] && continue

  # Acquire token
  while :; do
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
      0)  break ;;  # got a token
      1)
        sleep 1
        periodic_reset
        ;;
      99)
        echo "Lock starvation detected" >> "${LOG_DIR}/lock_errors.log"
        sleep 2
        ;;
    esac
  done

  if [[ "$DYNAMIC_RESOURCES" == "1" ]]; then
    read SLURM_TIME SLURM_MEM CPUS < <(estimate_resources "$ITEM")
  fi
  submit_job "$ITEM"


  if (( ++CLEANUP_COUNTER >= 20 )); then
    cleanup_successful_jobs &
    CLEANUP_COUNTER=0
  fi
done < "$ITEMS_FILE"

# Final cleanup
cleanup_successful_jobs

# Wait for all tokens to return => all jobs done
while [[ $(flock -x "$TOKEN_LOCK_FILE" -c "cat $TOKEN_FILE") -lt $MAX_JOBS ]]; do
  periodic_reset
  sleep 5
done

echo "🎉 All jobs completed!"
