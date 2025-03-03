#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --output=logs/kmer_jobs.log  # Single log file for all output
#SBATCH --error=logs/kmer_jobs.log   # Single log file for all errors
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Configuration (all in one file)
INPUT_DIR="/storage/home/mak6930/scratch/all/qc"       # Location of trimmed FASTQ files
OUTPUT_DIR="/storage/home/mak6930/scratch/all/kmer_counts"  # Where Jellyfish outputs go
CHECKPOINT_FILE="completed_kmer.txt"                  # File to track completed samples
HASH_SIZE="1G"                                       # Hash size parameter for Jellyfish
THREADS=4                                            # Number of threads for Jellyfish
KMER_RANGE="3 4 5 6 7 8"                             # K-mer sizes to process
BATCH_SIZE=50                                        # Number of jobs to submit at a time
MAX_RETRIES=3                                        # Max retries for Jellyfish failures
JOB_LOG="logs/kmer_jobs.log"                         # Single log file for all logs

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
        if [[ -f "$INPUT_DIR/${SAMPLE_NAME}_trimmed_1.fastq.gz" && -f "$INPUT_DIR/${SAMPLE_NAME}_trimmed_2.fastq.gz" ]]; then
            jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C -o "$JF_OUTPUT" \
                <(zcat "$INPUT_DIR/${SAMPLE_NAME}_trimmed_1.fastq.gz" "$INPUT_DIR/${SAMPLE_NAME}_trimmed_2.fastq.gz")
        elif [[ -f "$INPUT_DIR/${SAMPLE_NAME}_trimmed.fastq.gz" ]]; then
            jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C -o "$JF_OUTPUT" \
                <(zcat "$INPUT_DIR/${SAMPLE_NAME}_trimmed.fastq.gz")
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') No valid trimmed FASTQ found for sample $SAMPLE_NAME at k-mer $KMER_SIZE. Skipping." >> "$JOB_LOG"
            return 1
        fi

        # Check if Jellyfish succeeded
        if [[ $? -eq 0 ]]; then
            jellyfish dump -c -t -o "$TXT_OUTPUT" "$JF_OUTPUT"
            echo "$(date '+%Y-%m-%d %H:%M:%S') Finished processing sample $SAMPLE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"
            # Remove temporary Jellyfish binary file to save space
            rm -f "$JF_OUTPUT"
            return 0
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
    echo "$(date '+%Y-%m-%d %H:%M:%S') Starting processing for sample $SAMPLE_NAME" >> "$JOB_LOG"

    # Process each k-mer size in parallel using GNU Parallel
    export -f process_kmer
    export INPUT_DIR OUTPUT_DIR THREADS HASH_SIZE MAX_RETRIES JOB_LOG
    echo "$KMER_RANGE" | tr ' ' '\n' | parallel -j $THREADS "process_kmer '$SAMPLE_NAME' {}"

    # Mark the sample as completed if all k-mer sizes were processed successfully
    if [[ $? -eq 0 ]]; then
        echo "$SAMPLE_NAME" >> "$CHECKPOINT_FILE"
        echo "$(date '+%Y-%m-%d %H:%M:%S') Finished processing sample $SAMPLE_NAME" >> "$JOB_LOG"
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') Some k-mer sizes failed for sample $SAMPLE_NAME. Check logs for details." >> "$JOB_LOG"
    fi
}

# Get list of unprocessed samples
FILES_TO_PROCESS=()
for FILE in "$INPUT_DIR"/*_trimmed.fastq.gz "$INPUT_DIR"/*_trimmed_1.fastq.gz; do
    SAMPLE_NAME=$(basename "$FILE" | sed -E 's/(_trimmed_1|_trimmed_2|_trimmed)?\.fastq\.gz//')
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

# Submit jobs as a job array
echo "$(date '+%Y-%m-%d %H:%M:%S') Submitting jobs as a SLURM array..." >> "$JOB_LOG"
sbatch --array=1-${#FILES_TO_PROCESS[@]}%$BATCH_SIZE --job-name="kmer_array" --output="$JOB_LOG" --error="$JOB_LOG" \
    --time=02:00:00 --mem=8G --cpus-per-task=$THREADS --ntasks=1 \
    --wrap="$(declare -f process_sample process_kmer); export INPUT_DIR OUTPUT_DIR THREADS HASH_SIZE MAX_RETRIES KMER_RANGE JOB_LOG; process_sample '${FILES_TO_PROCESS[$SLURM_ARRAY_TASK_ID - 1]}'"

echo "$(date '+%Y-%m-%d %H:%M:%S') All jobs submitted!" >> "$JOB_LOG"
