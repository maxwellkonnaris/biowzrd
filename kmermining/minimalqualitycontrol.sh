#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

# Maximum number of simultaneous jobs
MAX_JOBS=50
CHECKPOINT_FILE="completed_qc.txt"

# Directories
INPUT_DIR="/storage/home/mak6930/scratch/all/fastq_data"
OUTPUT_DIR="/storage/home/mak6930/scratch/all/qc"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"

# Get list of unprocessed FASTQ samples
FILES_TO_PROCESS=()
for FILE in "$INPUT_DIR"/*.fastq.gz; do
    SAMPLE_NAME=$(basename "$FILE" | sed -E 's/(_1|_2)?\.fastq\.gz//')
    if ! grep -q "^$SAMPLE_NAME$" "$CHECKPOINT_FILE"; then
        FILES_TO_PROCESS+=("$SAMPLE_NAME")
    fi
done

TOTAL_FILES=${#FILES_TO_PROCESS[@]}
echo "Found $TOTAL_FILES unprocessed samples."

if [[ $TOTAL_FILES -eq 0 ]]; then
    echo "All samples have already been processed. Exiting."
    exit 0
fi

INDEX=0  # Track the current sample index

while [[ $INDEX -lt $TOTAL_FILES ]]; do
    # Check how many jobs are currently running
    RUNNING_JOBS=$(squeue -u $USER --name=fastq_* --format=%A | wc -l)
    
    # Submit new jobs if we have slots available
    while [[ $RUNNING_JOBS -lt $MAX_JOBS && $INDEX -lt $TOTAL_FILES ]]; do
        SAMPLE="${FILES_TO_PROCESS[$INDEX]}"
        JOB_SCRIPT="logs/job_${SAMPLE}.sh"
        
        # Generate SLURM script dynamically
        cat <<EOT > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=fastq_${SAMPLE}
#SBATCH --output=logs/${SAMPLE}.out
#SBATCH --error=logs/${SAMPLE}.err
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# Directories
INPUT_DIR="/storage/home/mak6930/scratch/all/fastq_data"
OUTPUT_DIR="/storage/home/mak6930/scratch/all/qc"
CHECKPOINT_FILE="completed_qc.txt"

# Parameters
THREADS=4

# Sample name
SAMPLE_NAME="$SAMPLE"
R1="\$INPUT_DIR/\${SAMPLE_NAME}_1.fastq.gz"
R2="\$INPUT_DIR/\${SAMPLE_NAME}_2.fastq.gz"

if [[ -f "\$R1" && -f "\$R2" ]]; then
    echo "Processing PAIRED-END sample: \$SAMPLE_NAME"

    fastp \
        -i "\$R1" \
        -I "\$R2" \
        -o "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_1.fastq.gz" \
        -O "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_2.fastq.gz" \
        --merged_out "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed.fastq.gz" \
        --html "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.html" \
        --json "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.json" \
        --thread "\$THREADS" \
        --detect_adapter_for_pe \
        --length_required 50 \
        --average_qual 20 \
        --cut_tail \
        --cut_tail_mean_quality 20 \
        --cut_tail_window_size 4 \
        --merge

else
    echo "Processing SINGLE-END sample: \$SAMPLE_NAME"

    fastp \
        -i "\$INPUT_DIR/\${SAMPLE_NAME}.fastq.gz" \
        -o "\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed.fastq.gz" \
        --html "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.html" \
        --json "\$OUTPUT_DIR/\${SAMPLE_NAME}_report.json" \
        --thread "\$THREADS" \
        --length_required 50 \
        --average_qual 20 \
        --cut_tail \
        --cut_tail_mean_quality 20 \
        --cut_tail_window_size 4

fi

# Mark this sample as completed
echo "\$SAMPLE_NAME" >> "\$CHECKPOINT_FILE"

echo "Finished processing \$SAMPLE_NAME!"
EOT

        # Submit the job
        sbatch "$JOB_SCRIPT"
        echo "Submitted job for sample: $SAMPLE"

        # Move to the next sample
        ((INDEX++))

        # Recalculate running jobs
        RUNNING_JOBS=$(squeue -u $USER --name=fastq_* --format=%A | wc -l)
    done

    echo "Currently running $RUNNING_JOBS jobs. Waiting before checking again..."

    # Wait before checking job count again
    sleep 60
done

echo "All jobs submitted!"

