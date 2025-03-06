#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --output=logs/kmer_jobs.log  # Single log file for all output
#SBATCH --error=logs/kmer_jobs.log   # Single log file for all errors
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Configuration 
INPUT_DIR="/storage/home/mak6930/scratch/all/qc"       # Location of FASTQ files
OUTPUT_DIR="/storage/home/mak6930/scratch/all/kmer_counts"  # Where Jellyfish outputs go
CHECKPOINT_FILE="completed_kmer.txt"                  # File to track completed k-mer sizes
HASH_SIZE="1G"                                       # Hash size parameter for Jellyfish
THREADS=4                                            # Number of threads for Jellyfish
KMER_RANGE="3 4 5 6 7 8"                             # K-mer sizes to process
MAX_JOBS=5                                           # Maximum number of jobs to run concurrently
MAX_RETRIES=3                                        # Max retries for Jellyfish failures
JOB_LOG="logs/kmer_jobs.log"                         # Single log file for all logs
ACTIVE_JOBS_FILE="logs/active_jobs.txt"              # Track jobs from this script only

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"
touch "$ACTIVE_JOBS_FILE"

# Get list of all FASTQ files in the input directory
FASTQ_FILES=("$INPUT_DIR"/*.fastq.gz)

# Filter out already processed k-mer sizes for each file
FILES_TO_PROCESS=()
for FILE in "${FASTQ_FILES[@]}"; do
    BASE_NAME=$(basename "$FILE")  # Keep full file name
    for KMER_SIZE in $KMER_RANGE; do
        if ! grep -q "^$BASE_NAME:$KMER_SIZE$" "$CHECKPOINT_FILE"; then
            FILES_TO_PROCESS+=("$BASE_NAME:$KMER_SIZE")
        fi
    done
done

TOTAL_FILES=${#FILES_TO_PROCESS[@]}
echo "$(date '+%Y-%m-%d %H:%M:%S') Found $TOTAL_FILES unprocessed k-mer sizes." >> "$JOB_LOG"

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') All k-mer sizes have been processed. Exiting." >> "$JOB_LOG"
    exit 0
fi

# Submit jobs, maintaining a maximum of MAX_JOBS running at any time
INDEX=0
while [[ $INDEX -lt $TOTAL_FILES ]]; do
    # Remove completed jobs from tracking file
    if [[ -s "$ACTIVE_JOBS_FILE" ]]; then
        while read -r job_id; do
            if ! squeue -j "$job_id" &>/dev/null; then
                sed -i "/^$job_id$/d" "$ACTIVE_JOBS_FILE"  # Remove completed job

                # Append job logs to the main log file and delete individual logs
                FILE_NAME=$(grep "$job_id" "$ACTIVE_JOBS_FILE" | awk '{print $2}')
                KMER_SIZE=$(grep "$job_id" "$ACTIVE_JOBS_FILE" | awk '{print $3}')
                if [[ -n "$FILE_NAME" && -n "$KMER_SIZE" ]]; then
                    echo "$(date '+%Y-%m-%d %H:%M:%S') Appending logs for job $job_id (file $FILE_NAME, k-mer $KMER_SIZE) to $JOB_LOG" >> "$JOB_LOG"
                    cat "logs/${FILE_NAME}_${KMER_SIZE}_kmer.out" >> "$JOB_LOG"
                    cat "logs/${FILE_NAME}_${KMER_SIZE}_kmer.err" >> "$JOB_LOG"
                    rm -f "logs/${FILE_NAME}_${KMER_SIZE}_kmer.out" "logs/${FILE_NAME}_${KMER_SIZE}_kmer.err"
                    echo "$(date '+%Y-%m-%d %H:%M:%S') Deleted individual logs for job $job_id (file $FILE_NAME, k-mer $KMER_SIZE)" >> "$JOB_LOG"
                fi

                # Delete the job script
                JOB_SCRIPT="logs/job_kmer_${FILE_NAME}_${KMER_SIZE}.sh"
                if [[ -f "$JOB_SCRIPT" ]]; then
                    rm -f "$JOB_SCRIPT"
                    echo "$(date '+%Y-%m-%d %H:%M:%S') Deleted job script: $JOB_SCRIPT" >> "$JOB_LOG"
                fi
            fi
        done < "$ACTIVE_JOBS_FILE"
    fi

    # Count running jobs from this script
    RUNNING_JOBS=$(wc -l < "$ACTIVE_JOBS_FILE")

    if [[ $RUNNING_JOBS -lt $MAX_JOBS ]]; then
        FILE_KMER="${FILES_TO_PROCESS[$INDEX]}"
        FILE_NAME=$(echo "$FILE_KMER" | cut -d':' -f1)
        KMER_SIZE=$(echo "$FILE_KMER" | cut -d':' -f2)
        JOB_SCRIPT="logs/job_kmer_${FILE_NAME}_${KMER_SIZE}.sh"

        # Create a job script for this file and k-mer size
       cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=kmer_${FILE_NAME}_${KMER_SIZE}
#SBATCH --output=logs/${FILE_NAME}_${KMER_SIZE}_kmer.out
#SBATCH --error=logs/${FILE_NAME}_${KMER_SIZE}_kmer.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

echo "$(date '+%Y-%m-%d %H:%M:%S') Starting processing for file $FILE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Define output files
JF_OUTPUT="$OUTPUT_DIR/${FILE_NAME}_kmer_${KMER_SIZE}.jf"
TXT_OUTPUT="$OUTPUT_DIR/${FILE_NAME}_kmer_counts_${KMER_SIZE}.txt"

# Debug: Print output file paths
echo "Output file (Jellyfish): $JF_OUTPUT" >> "$JOB_LOG"
echo "Output file (Text): $TXT_OUTPUT" >> "$JOB_LOG"

# Retry loop for Jellyfish
retry_count=0
while [[ $retry_count -lt $MAX_RETRIES ]]; do
    if [[ -f "$INPUT_DIR/$FILE_NAME" ]]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') Running Jellyfish for $FILE_NAME with k-mer size $KMER_SIZE (attempt $((retry_count + 1))/$MAX_RETRIES)" >> "$JOB_LOG"
        jellyfish count -m "$KMER_SIZE" -s "$HASH_SIZE" -t "$THREADS" -C -o "$JF_OUTPUT" \
            <(zcat "$INPUT_DIR/$FILE_NAME")
        
        # Check if Jellyfish succeeded
        if [[ $? -eq 0 ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') Jellyfish completed successfully for $FILE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"
            
            # Dump the k-mer counts into a human-readable format
            jellyfish dump -c -t -o "$TXT_OUTPUT" "$JF_OUTPUT"
            if [[ $? -eq 0 ]]; then
                echo "$(date '+%Y-%m-%d %H:%M:%S') K-mer counts dumped to $TXT_OUTPUT" >> "$JOB_LOG"
                
                # Append to checkpoint file only for this k-mer size
                echo "$FILE_NAME:$KMER_SIZE" >> "$CHECKPOINT_FILE"
                echo "$(date '+%Y-%m-%d %H:%M:%S') Added $FILE_NAME:$KMER_SIZE to checkpoint file." >> "$JOB_LOG"
                
                # Remove the temporary Jellyfish binary file to save space
                rm -f "$JF_OUTPUT"
                echo "$(date '+%Y-%m-%d %H:%M:%S') Removed temporary file: $JF_OUTPUT" >> "$JOB_LOG"
                break
            else
                echo "$(date '+%Y-%m-%d %H:%M:%S') Failed to dump k-mer counts for $FILE_NAME at k-mer $KMER_SIZE. Skipping." >> "$JOB_LOG"
                retry_count=$((retry_count + 1))
            fi
        else
            retry_count=$((retry_count + 1))
            echo "$(date '+%Y-%m-%d %H:%M:%S') Jellyfish failed for $FILE_NAME at k-mer $KMER_SIZE (attempt $retry_count/$MAX_RETRIES). Retrying..." >> "$JOB_LOG"
        fi
    else
        echo "$(date '+%Y-%m-%d %H:%M:%S') No valid FASTQ file found for $FILE_NAME at k-mer $KMER_SIZE. Skipping." >> "$JOB_LOG"
        break
    fi
done

echo "$(date '+%Y-%m-%d %H:%M:%S') Finished processing file $FILE_NAME with k-mer size $KMER_SIZE" >> "$JOB_LOG"
EOT

        # Make the job script executable
        chmod +x "$JOB_SCRIPT"

        # Submit the job and store the job ID
        JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $4}')
        echo "$JOB_ID $FILE_NAME $KMER_SIZE" >> "$ACTIVE_JOBS_FILE"

        INDEX=$((INDEX + 1))
        echo "$(date '+%Y-%m-%d %H:%M:%S') Submitted job $JOB_ID for file $FILE_NAME with k-mer size $KMER_SIZE. Running jobs: $RUNNING_JOBS" >> "$JOB_LOG"
    else
        sleep 10  # Wait before checking again
    fi
done

echo "$(date '+%Y-%m-%d %H:%M:%S') All jobs submitted!" >> "$JOB_LOG"
