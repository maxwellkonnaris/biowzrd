#!/bin/bash
#SBATCH --job-name=main_DL
#SBATCH --output=slurm_%A.out
#SBATCH --error=slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

#######################################
# Configuration
#######################################
MAX_JOBS=20  # Adjusted based on cluster capacity
WORKDIR="$(pwd)"
CHECKPOINT_FILE="${WORKDIR}/completed_accessions.txt"
ACCESSIONS_FILE="${WORKDIR}/run_accessions.txt"
COMBINED_METADATA="${WORKDIR}/metadata/combined_metadata.tsv"
DEBUG_LOG="${WORKDIR}/logs/debug.log"
TOKEN_FILE="${WORKDIR}/.job_tokens"  # Semaphore for job control
LOCK_FILE="${WORKDIR}/checkpoint.lock"
DEBUG_LOCK="${WORKDIR}/logs/debug.lock"

# Create necessary directories
mkdir -p "${WORKDIR}/jobs" "${WORKDIR}/logs" "${WORKDIR}/fastq_data" "${WORKDIR}/metadata"

# Initialize semaphore system
if [[ ! -f "$TOKEN_FILE" ]]; then
    echo $MAX_JOBS > "$TOKEN_FILE"
fi

# Ensure checkpoint file exists
touch "${CHECKPOINT_FILE}"

#######################################
# Function to Cleanup Successful Jobs
#######################################
cleanup_successful_jobs() {
    flock -x "$LOCK_FILE"  # Ensure atomic access to CHECKPOINT_FILE
    while read -r ACCESSION; do
        ACCESSION=$(echo "$ACCESSION" | xargs)  # Trim whitespace
        if [[ -z "$ACCESSION" ]]; then
            continue
        fi

        # Delete job script
        JOB_SCRIPT="${WORKDIR}/jobs/download_${ACCESSION}.sh"
        if [[ -f "$JOB_SCRIPT" ]]; then
            rm -f "$JOB_SCRIPT" && echo "✅ Deleted job script for $ACCESSION"
        fi

        # Delete logs
        LOG_OUT="${WORKDIR}/logs/${ACCESSION}.out"
        LOG_ERR="${WORKDIR}/logs/${ACCESSION}.err"
        if [[ -f "$LOG_OUT" ]]; then
            rm -f "$LOG_OUT" && echo "✅ Deleted log file: $LOG_OUT"
        fi
        if [[ -f "$LOG_ERR" ]]; then
            rm -f "$LOG_ERR" && echo "✅ Deleted log file: $LOG_ERR"
        fi
    done < "$CHECKPOINT_FILE"
    flock -u "$LOCK_FILE"
}

#######################################
# Function to Submit a Job for an Accession
#######################################
submit_job() {
    local ACCESSION=$1
    local JOB_SCRIPT="${WORKDIR}/jobs/download_${ACCESSION}.sh"

    # Trim whitespace and validate accession
    ACCESSION=$(echo "$ACCESSION" | xargs)
    [[ -z "$ACCESSION" ]] && return

    # Create the job script with token management
    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${ACCESSION}
#SBATCH --output=${WORKDIR}/logs/${ACCESSION}.out
#SBATCH --error=${WORKDIR}/logs/${ACCESSION}.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# Trap to release token on any exit
cleanup() {
    flock -x 200
    tokens=\$(< "$TOKEN_FILE")
    echo \$((tokens + 1)) > "$TOKEN_FILE"
    flock -u 200
}
exec 200>"${TOKEN_FILE}.lock"
trap cleanup EXIT

# Main processing logic (keep your existing workflow here)
export NCBI_API_KEY="9c9e61f98934800c1aab47c4066f394cde08"
WORKDIR="$WORKDIR"
CHECKPOINT_FILE="$CHECKPOINT_FILE"
LOCK_FILE="$LOCK_FILE"
SRA_FILES="${WORKDIR}/fastq_data/${ACCESSION}.sra"
FASTQ_DIR="${WORKDIR}/fastq_data"
METADATA_DIR="${WORKDIR}/metadata"
COMBINED_METADATA="$COMBINED_METADATA"
ACCESSION="$ACCESSION"
DEBUG_LOCK="$DEBUG_LOCK"



# Determine the provider based on the accession prefix
ACCESSION_PREFIX=\${ACCESSION:0:3}
echo "DEBUG: Accession = \$ACCESSION, Prefix = \$ACCESSION_PREFIX" >> "\${WORKDIR}/logs/\${ACCESSION}.out"

if [[ "\$ACCESSION_PREFIX" =~ ^(SRR|SRX|SRS|SRP)$ ]]; then
    PROVIDER="sra"
elif [[ "\$ACCESSION_PREFIX" =~ ^(ERR|ERX|ERS|ERP|DRR|DRX|DRS|DRP)$ ]]; then
    PROVIDER="ena"
else
    flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: Unknown accession type: \$ACCESSION' >> '\${WORKDIR}/logs/debug.log'";
    exit 1
fi

echo "DEBUG: Using provider: \$PROVIDER for accession: \$ACCESSION" >> "\${WORKDIR}/logs/\${ACCESSION}.out"

if [[ "\$PROVIDER" == "sra" ]]; then
    # Step 1: Prefetch
    echo "🔹 Prefetching SRA file for \$ACCESSION"
    prefetch "\$ACCESSION" --output-file "\$SRA_FILES" || { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} prefetch failed' >> '\${WORKDIR}/logs/debug.log'"; exit 1;}
    sleep 3  # Throttle after prefetch

    # Step 2: Convert to FASTQ
    echo "🔹 Converting SRA to FASTQ"
    fasterq-dump "\$SRA_FILES" \
        --outdir "\$FASTQ_DIR" \
        --threads 4 \
        --mem 8G \
        --split-3 || { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} fasterq-dump failed' >> '\${WORKDIR}/logs/debug.log'"; exit 1; }
    sleep 3  # Throttle after fasterq-dump

    # Step 3: Fetch Metadata with retries
    MAX_RETRIES=3
    RETRY_DELAY=3
    echo "🔹 Fetching metadata (max \$MAX_RETRIES attempts)"
    for ((i=1; i<=\$MAX_RETRIES; i++)); do
        echo "Attempt \$i/3..."
        esearch -db sra -query "\"$ACCESSION\"" | efetch -db sra -format runinfo > "\${METADATA_DIR}/\${ACCESSION}-run-info.csv"
        if [[ \$? -eq 0 && -s "\${METADATA_DIR}/\${ACCESSION}-run-info.csv" ]]; then
        echo "Metadata fetched successfully"
        # Convert CSV to TSV and merge with combined metadata
            (
                flock -x 200  # Use file lock for combined metadata
                TSV_FILE="\${METADATA_DIR}/\${ACCESSION}-run-info.tsv"
                
                # Convert CSV to TSV
                tr ',' '\t' < "\${METADATA_DIR}/\${ACCESSION}-run-info.csv" > "\$TSV_FILE"
                
                # Append to combined metadata
                if [[ ! -f "\$COMBINED_METADATA" ]]; then
                    head -n 1 "\$TSV_FILE" > "\$COMBINED_METADATA"
                fi
                tail -n +2 "\$TSV_FILE" >> "\$COMBINED_METADATA"
                
                # Cleanup
                rm "\${METADATA_DIR}/\${ACCESSION}-run-info.csv"
                echo "🔹 Metadata appended to combined file"
            ) 200>"\$LOCK_FILE"
            
            break
        else
            flock -x "\$DEBUG_LOCK" bash -c "echo '⚠️ WARNING: \${ACCESSION} Metadata attempt \$i failed' >> '\${WORKDIR}/logs/debug.log'"
            sleep \$((RETRY_DELAY * i))
        fi
    done

    # Final check for metadata
    if [[ ! -s "\${METADATA_DIR}/\${ACCESSION}-run-info.csv" ]]; then
        flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} All metadata attempts failed' >> '\${WORKDIR}/logs/debug.log'"
    fi

elif [[ "\$PROVIDER" == "ena" ]]; then
    # Use fastq-dl for ENA accessions (includes metadata)
    echo "🔹 Downloading FASTQ and metadata for \$ACCESSION"
    fastq-dl --accession "\$ACCESSION" \
             --provider "\$PROVIDER" \
             --cpus 4 \
             --prefix "\$ACCESSION" \
             --outdir "\$FASTQ_DIR" || { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} fastq-dl failed' >> '\${WORKDIR}/logs/debug.log'"; exit 1; }

    # Process ENA metadata
    ENA_METADATA="\${FASTQ_DIR}/\${ACCESSION}-run-info.tsv"
    if [[ -f "\$ENA_METADATA" ]]; then
        # Use file locking for atomic operations
        (
            flock -x 200
            # Create combined file with header if missing
            if [[ ! -f "\$COMBINED_METADATA" ]]; then
                if ! head -n 1 "\$ENA_METADATA" > "\$COMBINED_METADATA"; then
                    echo "❌ ERROR: Failed to create combined metadata header" >&2
                fi
            fi
            
            # Append data rows
            if ! tail -n +2 "\$ENA_METADATA" >> "\$COMBINED_METADATA"; then
                echo "❌ ERROR: Failed to append metadata for \$ACCESSION" >&2
            fi
            
            # Move metadata file with cleanup
            if ! mv "\$ENA_METADATA" "\${METADATA_DIR}/"; then
                echo "❌ ERROR: Failed to move metadata file for \$ACCESSION" >&2
                # Cleanup partial data from combined file
                if [[ -f "\$COMBINED_METADATA" ]]; then
                    grep -v "^\${ACCESSION}" "\$COMBINED_METADATA" > "\$COMBINED_METADATA.tmp" && \
                    mv "\$COMBINED_METADATA.tmp" "\$COMBINED_METADATA"
                fi
            fi
            
            echo "🔹 Metadata processed and moved"
            
        ) 200>"\$LOCK_FILE"
    else
        # Log warning with timestamp and clean up any potential partial files
        (
            flock -x "\$DEBUG_LOCK"
            echo "\$(date '+%F %T') ⚠️ WARNING: \${ACCESSION} ENA metadata file not found" >> "\${WORKDIR}/logs/debug.log"
            # Remove potential empty file artifacts
            rm -f "\${METADATA_DIR}/\${ACCESSION}-run-info.tsv"
        )
    fi

else
    flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} Unsupported provider: \$PROVIDER' >> '\${WORKDIR}/logs/debug.log'";
    exit 1  # Important: Exit on unknown provider
fi

# Verify FASTQ files and compress
shopt -s nullglob
for FILE in "\${FASTQ_DIR}"/"\${ACCESSION}"*.fastq; do
    if [[ -f "\$FILE" && "\$FILE" != *.gz ]]; then
        echo "🔹 Gzipping: \$FILE"
        gzip "\$FILE" || { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} Gzip failed' >> '\${WORKDIR}/logs/debug.log'"; }
    fi
done

# Check for gzipped FASTQs
ALL_GZ_FASTQS=("\${FASTQ_DIR}"/"\${ACCESSION}"*.fastq.gz)
[[ \${#ALL_GZ_FASTQS[@]} -eq 0 ]] && { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} No FASTQ files' >> '\${WORKDIR}/logs/debug.log'"; exit 1; }


# Record success in checkpoint
flock -x "$LOCK_FILE" bash -c "echo \$ACCESSION >> \$CHECKPOINT_FILE"

# Cleanup SRA files
if ls "\$SRA_FILES" 1>/dev/null 2>&1; then
    rm -f "\$SRA_FILES" || { flock -x "\$DEBUG_LOCK" bash -c "echo '❌ ERROR: \${ACCESSION} Failed to delete SRA files' >> '\${WORKDIR}/logs/debug.log'"; exit 1; }
fi

echo "✅ Successfully processed \$ACCESSION"
EOF

    chmod +x "$JOB_SCRIPT"
    sbatch "$JOB_SCRIPT" && echo "✅ Submitted $ACCESSION"
}

#######################################
# Main Submission Loop with Token Control
#######################################
while read -r ACCESSION; do
    # Skip processed accessions
    grep -Fxq "$ACCESSION" "$CHECKPOINT_FILE" && continue

    # Acquire token with non-blocking flock
    while :; do
        flock -x "${TOKEN_FILE}.lock" -c "
            tokens=\$(< "$TOKEN_FILE")
            if (( tokens > 0 )); then
                echo \$((tokens - 1)) > "$TOKEN_FILE"
                exit 0
            fi
            exit 1
        " && break
        sleep $(( RANDOM % 5 + 1 ))
    done

    submit_job "$ACCESSION"
    sleep 0.5  # Throttle job submissions

    # Periodically cleanup successful jobs
    if (( RANDOM % 10 == 0 )); then  # Cleanup ~10% of the time
        cleanup_successful_jobs
    fi
done < "$ACCESSIONS_FILE"

# Final cleanup after all jobs are submitted
cleanup_successful_jobs

# Wait for final jobs to complete
while [[ $(flock -x "${TOKEN_FILE}.lock" -c "cat $TOKEN_FILE") -lt $MAX_JOBS ]]; do
    sleep 10
done

echo "🎉 All jobs completed!"
