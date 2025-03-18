#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=fastq_main
#SBATCH --output=logs/slurm_main_%A.out
#SBATCH --error=logs/slurm_main_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

#######################################
# Configuration
#######################################
MAX_JOBS=50  # Maximum number of concurrent jobs
WORKDIR="$(pwd)"
ACCESSIONS_FILE="${WORKDIR}/run_accessions.txt"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
TOKEN_FILE="${WORKDIR}/.job_tokens"
TOKEN_LOCK_FILE="${WORKDIR}/.job_tokens.lock"
LOG_DIR="${WORKDIR}/logs"
FASTQ_DIR="${WORKDIR}/fastq_files"
SCRIPT_DIR="${WORKDIR}/jobs"
PYTHON_SCRIPT="${WORKDIR}/download_fastq.py"

# Create necessary directories
mkdir -p "$LOG_DIR" "$FASTQ_DIR" "$SCRIPT_DIR"

# Ensure checkpoint file exists
touch "$CHECKPOINT_FILE"

# Initialize token system if missing
if [[ ! -f "$TOKEN_FILE" ]]; then
    echo "$MAX_JOBS" > "$TOKEN_FILE"
fi

#######################################
# Function: Submit a Job for an Accession
#######################################
submit_job() {
    local ACCESSION=$1
    local JOB_SCRIPT="${SCRIPT_DIR}/download_${ACCESSION}.sh"

    # Trim whitespace and validate accession
    ACCESSION=$(echo "$ACCESSION" | xargs)
    [[ -z "$ACCESSION" ]] && return

    # Create Slurm job script
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=${LOG_DIR}/${ACCESSION}.out
#SBATCH --error=${LOG_DIR}/${ACCESSION}.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

set -euo pipefail

# Run the Python script for this accession
python3 "$PYTHON_SCRIPT" "$ACCESSION"

# Mark as completed
echo "$ACCESSION" >> "$CHECKPOINT_FILE"

# Release token
(
    flock -x 9
    current_tokens=\$(< "$TOKEN_FILE")
    new_tokens=\$((current_tokens + 1))
    echo "\$new_tokens" > "$TOKEN_FILE"
) 9>"$TOKEN_LOCK_FILE"
EOF

    chmod +x "$JOB_SCRIPT"
    sbatch "$JOB_SCRIPT" && echo "✅ Submitted $ACCESSION"
}

#######################################
# Main Job Submission Loop
#######################################
while read -r ACCESSION; do
    # Skip already completed jobs
    grep -Fxq "$ACCESSION" "$CHECKPOINT_FILE" && continue

    # Acquire a token
    while :; do
        (
            flock -x 9 || exit 99  # Ensure atomic access
            current_tokens=$(< "$TOKEN_FILE")
            if (( current_tokens > 0 )); then
                echo $((current_tokens - 1)) > "$TOKEN_FILE"
                exit 0
            fi
            exit 1
        ) 9>"$TOKEN_LOCK_FILE"

        case $? in
            0) break ;;  # Token acquired
            1) sleep 2 ;; # Wait and retry
            99) echo "⚠ Lock issue detected" >> "$LOG_DIR/lock_errors.log"; sleep 5 ;;
        esac
    done

    # Submit the job
    submit_job "$ACCESSION"

done < "$ACCESSIONS_FILE"

# Wait for final jobs to finish
while [[ $(< "$TOKEN_FILE") -lt "$MAX_JOBS" ]]; do
    sleep 5
done

echo "🎉 All jobs completed!"
