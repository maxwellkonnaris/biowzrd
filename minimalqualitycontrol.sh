#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

###############################
#    USER DEFINED PARAMETERS  #
###############################

# Max concurrency (simultaneous jobs)
MAX_JOBS=10

CHECKPOINT_FILE="completed_qc.txt"
INPUT_DIR="/storage/home/mak6930/scratch/all/fastq_data"
OUTPUT_DIR="/storage/home/mak6930/scratch/all/qc"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"

###############################
#  COLLECT UNPROCESSED SAMPLES
###############################
FILES_TO_PROCESS=()

# Loop over all fastq.gz files, handling both single- and paired-end
for FILE in "$INPUT_DIR"/*.fastq.gz; do
    # e.g. "SAMPLE_1.fastq.gz" or "SAMPLE_2.fastq.gz" or "SAMPLE.fastq.gz"
    BASENAME=$(basename "$FILE" .fastq.gz)  # e.g. "SAMPLE_1" or "SAMPLE_2" or "SAMPLE"

    # If BASENAME ends with "_1" or "_2", strip that off
    if [[ $BASENAME =~ _[12]$ ]]; then
        SAMPLE_NAME="${BASENAME%_*}"  # remove last underscore and digit
    else
        # Single-end, or no trailing "_1"/"_2"
        SAMPLE_NAME="$BASENAME"
    fi

    # Skip if it's already in the checkpoint
    if grep -q "^$SAMPLE_NAME$" "$CHECKPOINT_FILE"; then
        continue
    fi

    # Avoid duplicates: if we've already queued this sample, skip
    if [[ " ${FILES_TO_PROCESS[@]} " =~ " $SAMPLE_NAME " ]]; then
        continue
    fi

    FILES_TO_PROCESS+=( "$SAMPLE_NAME" )
done

TOTAL_FILES=${#FILES_TO_PROCESS[@]}
echo "Found $TOTAL_FILES unprocessed samples."

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "All samples have already been processed. Exiting."
    exit 0
fi

###############################
#  JOB SUBMISSION LOOP
###############################

INDEX=0
while [[ $INDEX -lt $TOTAL_FILES ]]; do

    SAMPLE="${FILES_TO_PROCESS[$INDEX]}"

    ###############################
    #  CHECK CURRENT JOB COUNT
    ###############################
    while true; do
        # Count how many 'fastq_' jobs are currently running
        RUNNING_JOBS=$(squeue -u "$USER" -o "%A %j" | grep -c "fastq_")
        
        if [[ $RUNNING_JOBS -lt $MAX_JOBS ]]; then
            # We have a free slot to submit
            break
        else
            echo "Currently running $RUNNING_JOBS jobs. Max is $MAX_JOBS. Waiting..."
            sleep 30
        fi
    done

    ###############################
    #  CREATE & SUBMIT JOB SCRIPT
    ###############################
    JOB_SCRIPT="logs/job_${SAMPLE}.sh"
    cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${SAMPLE}
#SBATCH --output=logs/${SAMPLE}.out
#SBATCH --error=logs/${SAMPLE}.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

# Directories
INPUT_DIR="$INPUT_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
CHECKPOINT_FILE="$CHECKPOINT_FILE"

# Threads for fastp
THREADS=8
SAMPLE_NAME="$SAMPLE"

# Potential paired-end files
R1="\$INPUT_DIR/\${SAMPLE_NAME}_1.fastq.gz"
R2="\$INPUT_DIR/\${SAMPLE_NAME}_2.fastq.gz"

if [[ -f "\$R1" && -f "\$R2" ]]; then
    echo "Processing PAIRED-END sample: \$SAMPLE_NAME"

    fastp \\
        -i "\$R1" \\
        -I "\$R2" \\
        -o "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_1.fastq.gz" \\
        -O "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_2.fastq.gz" \\
        --merged_out "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed.fastq.gz" \\
        --html "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.html" \\
        --json "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.json" \\
        --thread "\$THREADS" \\
        --detect_adapter_for_pe \\
        --length_required 50 \\
        --average_qual 20 \\
        --cut_tail \\
        --cut_tail_mean_quality 20 \\
        --cut_tail_window_size 4 \\
        --merge
elif [[ -f "\$INPUT_DIR/\${SAMPLE_NAME}.fastq.gz" ]]; then
    echo "Processing SINGLE-END sample: \$SAMPLE_NAME"

    fastp \\
        -i "\$INPUT_DIR/\${SAMPLE_NAME}.fastq.gz" \\
        -o "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed.fastq.gz" \\
        --html "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.html" \\
        --json "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.json" \\
        --thread "\$THREADS" \\
        --length_required 50 \\
        --average_qual 20 \\
        --cut_tail \\
        --cut_tail_mean_quality 20 \\
        --cut_tail_window_size 4
else
    echo "ERROR: Could not find appropriate files for sample: \$SAMPLE_NAME"
    echo "Skipping \$SAMPLE_NAME" 1>&2
    exit 1
fi

# Mark this sample as completed
echo "\$SAMPLE_NAME" >> "\$CHECKPOINT_FILE"
echo "Finished processing \$SAMPLE_NAME!"
EOT

    sbatch "$JOB_SCRIPT"
    echo "Submitted job for sample: $SAMPLE"

    # Move to the next sample
    ((INDEX++))
done

###############################
#  (OPTIONAL) WAIT FOR ALL JOBS
###############################
# If you want this script to exit only after all jobs finish, uncomment this:

# echo "All jobs submitted. Waiting for all 'fastq_' jobs to complete..."
# while true; do
#     RUNNING_JOBS=$(squeue -u "$USER" -o "%A %j" | grep -c "fastq_")
#     if [[ $RUNNING_JOBS -eq 0 ]]; then
#         echo "All 'fastq_' jobs have completed."
#         break
#     else
#         echo "Still \$RUNNING_JOBS jobs running... checking again in 30 seconds."
#         sleep 30
#     fi
# done

echo "Done submitting jobs."

