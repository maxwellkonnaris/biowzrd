#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --output=logs/kmer_jobs.log  # Single log file for all output
#SBATCH --error=logs/kmer_jobs.log   # Single log file for all errors
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Configuration (all in one file)
INPUT_DIR="/storage/home/mak6930/scratch/all/qc"       # Location of FASTQ files
OUTPUT_DIR="/storage/home/mak6930/scratch/all/kmer_counts"  # Where Jellyfish outputs go
CHECKPOINT_FILE="completed_kmer.txt"                  # File to track completed samples
HASH_SIZE="1G"                                       # Hash size parameter for Jellyfish
THREADS=4                                            # Number of threads for Jellyfish
KMER_RANGE="3 4 5 6 7 8"                             # K-mer sizes to process
BATCH_SIZE=50                                        # Number of jobs to submit at a time
MAX_RETRIES=3                                        # Max retries for Jellyfish failures
JOB_LOG="logs/kmer_jobs.log"                         # Single log file for all logs
MAX_JOBS=1000                                        # Maximum number of child jobs to submit

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"

# Function to process a single sample and k-mer size
process_kmer() {
    SAMPLE_NAME="$1"
    KMER_SIZE="$2"
    echo "$(date '+%Y-%m-%d %H:%M:%S') Processing sample $SAMPLE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"

    JF_OUTPUT="$OUTPUT_DIR/${SAMPLE_NAME}_kmer_${KMER_SIZE}.jf"
    TXT_OUTPUT="$OUTPUT_DIR/${SAMPLE_NAME}_kmer_counts_${KMER_SIZE}.txt"

    # Retry loop for Jellyfish
    retry_count=0
    while [[ $retry_count -lt $MAX_RETRIES ]]; do
        # Check for paired-end files (Sample_name_1.fastq.gz and Sample_name_2.fastq.gz)
        if [[ -f "$INPUT_DIR/${SAMPLE_NAME}_1.fastq.gz" && -f "$INPUT_DIR/${SAMPLE_NAME}_2.fastq.gz" ]]; then
            jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C -o "$JF_OUTPUT" \
                <(zcat "$INPUT_DIR/${SAMPLE_NAME}_1.fastq.gz" "$INPUT_DIR/${SAMPLE_NAME}_2.fastq.gz")
        # Check for single-end files (Sample_name.fastq.gz)
        elif [[ -f "$INPUT_DIR/${SAMPLE_NAME}.fastq.gz" ]]; then
            jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C -o "$JF_OUTPUT" \
                <(zcat "$INPUT_DIR/${SAMPLE_NAME}.fastq.gz")
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') No valid FASTQ files found for sample $SAMPLE_NAME at k-mer $KMER_SIZE. Skipping." >> "$JOB_LOG"
            return 1
        fi

        # Check if Jellyfish succeeded
        if [[ $? -eq 0 ]]; then
            # Dump the k-mer counts into a human-readable format
            jellyfish dump -c -t -o "$TXT_OUTPUT" "$JF_OUTPUT"
            if [[ $? -eq 0 ]]; then
                echo "$(date '+%Y-%m-%d %H:%M:%S') Finished processing sample $SAMPLE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"
                # Remove the temporary Jellyfish binary file to save space
                rm -f "$JF_OUTPUT"
                echo "$(date '+%Y-%m-%d %H:%M:%S') Removed temporary file: $JF_OUTPUT" >> "$JOB_LOG"
                return 0
            else
                echo "$(date '+%Y-%m-%d %H:%M:%S') Failed to dump k-mer counts for sample $SAMPLE_NAME at k-mer $KMER_SIZE. Skipping." >> "$JOB_LOG"
                return 1
            fi
        else
            retry_count=$((retry_count + 1))
            echo "$(date '+%Y-%m-%d %H:%M:%S') Jellyfish failed for sample $SAMPLE_NAME at k-mer $KMER_SIZE (attempt $retry_count/$MAX_RETRIES). Retrying..." >> "$JOB_LOG"
        fi
    done

    echo "$(date '+%Y-%m-%d %H:%M:%S') Jellyfish failed for sample $SAMPLE_NAME at k-mer $KMER_SIZE after $MAX_RETRIES attempts. Skipping." >> "$JOB_LOG"
    return 1
}

# Function to process all k-mer sizes for a single sample
process_sample() {
    SAMPLE_NAME="$1"
    JOB_SCRIPT="logs/job_kmer_${SAMPLE_NAME}.sh"
    echo "$(date '+%Y-%m-%d %H:%M:%S') Starting processing for sample $SAMPLE_NAME" >> "$JOB_LOG"

    # Create a job script for this sample
    cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=kmer_${SAMPLE_NAME}
#SBATCH --output=logs/${SAMPLE_NAME}_kmer.out
#SBATCH --error=logs/${SAMPLE_NAME}_kmer.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --dependency=afterok:$SLURM_JOB_ID  # Add dependency on the parent job

# Process each k-mer size in parallel using GNU Parallel
export -f process_kmer
export INPUT_DIR OUTPUT_DIR THREADS HASH_SIZE MAX_RETRIES JOB_LOG
echo "$KMER_RANGE" | tr ' ' '\n' | parallel -j $THREADS "process_kmer '$SAMPLE_NAME' {}"

# Mark the sample as completed if all k-mer sizes were processed successfully
if [[ \$? -eq 0 ]]; then
    echo "$SAMPLE_NAME" >> "$CHECKPOINT_FILE"
    echo "$(date '+%Y-%m-%d %H:%M:%S') Finished processing sample $SAMPLE_NAME" >> "$JOB_LOG"
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') Some k-mer sizes failed for sample $SAMPLE_NAME. Check logs for details." >> "$JOB_LOG"
fi

# Remove the job script after completion
rm -f "$JOB_SCRIPT"
echo "$(date '+%Y-%m-%d %H:%M:%S') Removed job script: $JOB_SCRIPT" >> "$JOB_LOG"
EOT

    # Make the job script executable
    chmod +x "$JOB_SCRIPT"

    # Submit the job script
    sbatch "$JOB_SCRIPT"
}

# Get list of all FASTQ files in the input directory
FASTQ_FILES=("$INPUT_DIR"/*.fastq.gz)

# Extract unique sample names
declare -A SAMPLE_NAMES
for FILE in "${FASTQ_FILES[@]}"; do
    # Extract base name (e.g., Sample_name_1.fastq.gz -> Sample_name)
    BASE_NAME=$(basename "$FILE" .fastq.gz)
    SAMPLE_NAME="${BASE_NAME%_[12]}"  # Remove _1 or _2 suffix for paired-end files
    SAMPLE_NAMES["$SAMPLE_NAME"]=1
done

# Convert associative array keys to a list of unique sample names
UNIQUE_SAMPLES=("${!SAMPLE_NAMES[@]}")

# Filter out already processed samples
FILES_TO_PROCESS=()
for SAMPLE_NAME in "${UNIQUE_SAMPLES[@]}"; do
    if ! grep -q "^$SAMPLE_NAME$" "$CHECKPOINT_FILE"; then
        FILES_TO_PROCESS+=("$SAMPLE_NAME")
    fi
done

TOTAL_FILES=${#FILES_TO_PROCESS[@]}
echo "$(date '+%Y-%m-%d %H:%M:%S') Found $TOTAL_FILES unprocessed samples." >> "$JOB_LOG"

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') All samples have been processed. Exiting." >> "$JOB_LOG"
    exit 0
fi

# Submit jobs in batches
JOB_COUNT=0
for (( i=0; i<TOTAL_FILES; i+=BATCH_SIZE )); do
    # Check if the parent job is still running
    if ! squeue -j $SLURM_JOB_ID &> /dev/null; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Parent job canceled. Exiting." >> "$JOB_LOG"
        exit 1
    fi

    # Check if the maximum number of jobs has been reached
    if [[ $JOB_COUNT -ge $MAX_JOBS ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Reached maximum number of jobs ($MAX_JOBS). Exiting." >> "$JOB_LOG"
        exit 0
    fi

    BATCH=("${FILES_TO_PROCESS[@]:i:BATCH_SIZE}")
    for SAMPLE in "${BATCH[@]}"; do
        process_sample "$SAMPLE"
        JOB_COUNT=$((JOB_COUNT + 1))
    done
    echo "$(date '+%Y-%m-%d %H:%M:%S') Submitted batch of ${#BATCH[@]} jobs." >> "$JOB_LOG"
    sleep 5  # Short pause to prevent overwhelming the scheduler
done

echo "$(date '+%Y-%m-%d %H:%M:%S') All jobs submitted!" >> "$JOB_LOG"
