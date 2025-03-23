#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=QC
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

#######################################
# Configuration
#######################################
MAX_JOBS=20
WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_qc.txt"
METADATA_FILE="${WORKDIR}/sample_types.txt"
INPUT_DIR="/storage/home/mak6930/scratch/all/fastq_data"
OUTPUT_DIR="/storage/home/mak6930/scratch/all/qc"
DEBUG_LOG="${WORKDIR}/logs/debug.log"
TOKEN_FILE="${WORKDIR}/.job_tokens"
CHECKPOINT_LOCK_FILE="${WORKDIR}/checkpoint.lock"
TOKEN_LOCK_FILE="${WORKDIR}/.job_tokens.lock"
PERIODIC_RESET_LOG="${WORKDIR}/logs/periodic_reset.log"
JOB_NAME_PREFIX="qc_"

# Create directories
mkdir -p "${OUTPUT_DIR}" "${WORKDIR}/logs" "${WORKDIR}/jobs"

# Initialize semaphore system
echo "$MAX_JOBS" > "$TOKEN_FILE"
touch "$CHECKPOINT_FILE"

#######################################
# Load Sample Type Metadata
#######################################
declare -A SAMPLE_TYPES
while IFS=$'\t' read -r sample type; do
    SAMPLE_TYPES["$sample"]="$type"
done < "$METADATA_FILE"

#######################################
# Periodic Token Reset (From Reference)
#######################################
periodic_reset() {
  {
    flock -x 9 || return 1
    current_tokens=$(cat "$TOKEN_FILE" 2>/dev/null || echo "$MAX_JOBS")
    
    # Count actual running jobs
    RUNNING_JOBS=$(squeue --noheader --format '%j' | grep -c "^${JOB_NAME_PREFIX}")
    tokens_out=$((MAX_JOBS - current_tokens))

    if (( tokens_out != RUNNING_JOBS )); then
      echo "[$(date)] Resetting tokens: had $current_tokens, running $RUNNING_JOBS" >> "$PERIODIC_RESET_LOG"
      new_tokens=$((MAX_JOBS - RUNNING_JOBS))
      echo "$new_tokens" > "$TOKEN_FILE"
    fi
  } 9<>"$TOKEN_LOCK_FILE"
}

#######################################
# Cleanup Successful Jobs (From Reference)
#######################################
cleanup_successful_jobs() {
  flock -x "$CHECKPOINT_LOCK_FILE" || return 1
  
  tail -n 100 "$CHECKPOINT_FILE" | while read -r SAMPLE; do
    # Delete job scripts
    find "${WORKDIR}/jobs" -name "job_${SAMPLE}.sh" -delete
    
    # Delete logs older than 7 days
    find "${WORKDIR}/logs" -name "${SAMPLE}.*" -mtime +7 -delete
  done

  flock -u "$CHECKPOINT_LOCK_FILE"
}

#######################################
# Main Submission Loop with Token Control
#######################################
CLEANUP_COUNTER=0
while read -r SAMPLE; do
    # Skip completed or invalid samples
    grep -Fxq "$SAMPLE" "$CHECKPOINT_FILE" && continue
    [[ -z "${SAMPLE_TYPES[$SAMPLE]}" ]] && continue

    # Token acquisition loop
    while :; do
        (
            flock -x 9 || exit 99
            current_tokens=$(< "$TOKEN_FILE")
            if (( current_tokens > 0 )); then
                echo $((current_tokens - 1)) > "$TOKEN_FILE"
                break
            fi
            exit 1
        ) 9<>"$TOKEN_LOCK_FILE"

        case $? in
            0) break;;
            1) sleep 5; periodic_reset;;
            99) sleep 10;;  # Lock contention
        esac
    done

    #######################################
    # Create Job Script with Token Release
    #######################################
    JOB_SCRIPT="${WORKDIR}/jobs/job_${SAMPLE}.sh"
    SAMPLE_TYPE="${SAMPLE_TYPES[$SAMPLE]}"
    
    cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=${JOB_NAME_PREFIX}${SAMPLE}
#SBATCH --output=logs/${SAMPLE}.out
#SBATCH --error=logs/${SAMPLE}.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

set -eo pipefail

# Token release mechanism
release_token() {
    (
        flock -x 9 || exit 1
        tokens=\$(< "$TOKEN_FILE")
        echo \$((tokens + 1)) > "$TOKEN_FILE"
    ) 9<>"$TOKEN_LOCK_FILE"
}
trap 'release_token' EXIT TERM INT

# Main processing
R1="${INPUT_DIR}/${SAMPLE}_1.fastq.gz"
R2="${INPUT_DIR}/${SAMPLE}_2.fastq.gz"

if [[ -f "\$R1" && -f "\$R2" ]]; then
    fastp \\
        -i "\$R1" -I "\$R2" \\
        -o "${OUTPUT_DIR}/${SAMPLE}_trimmed_1.fastq.gz" \\
        -O "${OUTPUT_DIR}/${SAMPLE}_trimmed_2.fastq.gz" \\
        $(if [[ "$SAMPLE_TYPE" == "16S" ]]; then
            echo "--trim_poly_g --trim_poly_x --cut_front --cut_tail"
            echo "--cut_window_size 4 --cut_mean_quality 20"
            echo "--qualified_quality_phred 20 --length_required 100 --n_base_limit 0"
        else
            echo "--trim_poly_g --trim_poly_x --cut_front --cut_tail"
            echo "--cut_window_size 4 --cut_mean_quality 15"
            echo "--qualified_quality_phred 15 --length_required 50 --n_base_limit 5"
        fi) \\
        --thread 8 \\
        --detect_adapter_for_pe \\
        --json "${OUTPUT_DIR}/${SAMPLE}_report.json" \\
        --html "${OUTPUT_DIR}/${SAMPLE}_report.html"
elif [[ -f "${INPUT_DIR}/${SAMPLE}.fastq.gz" ]]; then
    fastp \\
        -i "${INPUT_DIR}/${SAMPLE}.fastq.gz" \\
        -o "${OUTPUT_DIR}/${SAMPLE}_trimmed.fastq.gz" \\
        $(if [[ "$SAMPLE_TYPE" == "16S" ]]; then
            echo "--trim_poly_g --trim_poly_x --cut_front --cut_tail"
            echo "--cut_window_size 4 --cut_mean_quality 20"
            echo "--qualified_quality_phred 20 --length_required 100 --n_base_limit 0"
        else
            echo "--trim_poly_g --trim_poly_x --cut_front --cut_tail"
            echo "--cut_window_size 4 --cut_mean_quality 15"
            echo "--qualified_quality_phred 15 --length_required 50 --n_base_limit 5"
        fi) \\
        --thread 8 \\
        --json "${OUTPUT_DIR}/${SAMPLE}_report.json" \\
        --html "${OUTPUT_DIR}/${SAMPLE}_report.html"
else
    echo "Missing files for ${SAMPLE}" >&2
    exit 1
fi

# Record completion
flock -x "$CHECKPOINT_LOCK_FILE" -c "echo '$SAMPLE' >> '$CHECKPOINT_FILE'"
EOT

    sbatch "$JOB_SCRIPT"
    echo "Submitted $SAMPLE (${SAMPLE_TYPE})"

    # Periodic cleanup
    ((CLEANUP_COUNTER++))
    if (( CLEANUP_COUNTER >= 20 )); then
        cleanup_successful_jobs &
        CLEANUP_COUNTER=0
    fi
done < <(printf "%s\n" "${!SAMPLE_TYPES[@]}")

# Final cleanup and wait
cleanup_successful_jobs
while [[ $(flock -x "$TOKEN_LOCK_FILE" -c "cat $TOKEN_FILE") -lt $MAX_JOBS ]]; do
    sleep 5
    periodic_reset
done

echo "✅ All QC jobs completed successfully!"
