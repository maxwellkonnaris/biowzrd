#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

###############################
#    USER DEFINED PARAMETERS  #
###############################

# Max concurrency (simultaneous jobs)
MAX_JOBS=20

# Path to BBMap's built-in adapter file (verify this path!)
ADAPTER_FILE="$CONDA_PREFIX/opt/bbmap/resources/adapters.fa"

CHECKPOINT_FILE="completed_qc.txt"
INPUT_DIR="/storage/home/mak6930/scratch/all/fastq_data"
OUTPUT_DIR="/storage/home/mak6930/scratch/all/qc"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Create checkpoint file if it doesn’t exist
touch "$CHECKPOINT_FILE"

###############################
#  COLLECT UNPROCESSED SAMPLES
###############################
FILES_TO_PROCESS=()

# Loop over all fastq.gz files, handling both single- and paired-end
for FILE in "$INPUT_DIR"/*.fastq.gz; do
    BASENAME=$(basename "$FILE" .fastq.gz)

    # If BASENAME ends with "_1" or "_2", strip that off
    if [[ $BASENAME =~ _[12]$ ]]; then
        SAMPLE_NAME="${BASENAME%_*}"
    else
        SAMPLE_NAME="$BASENAME"
    fi

    # Skip if already processed
    if grep -q "^$SAMPLE_NAME$" "$CHECKPOINT_FILE"; then
        continue
    fi

    # Avoid duplicates
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
#  FUNCTION: DETERMINE MEMORY
###############################
get_memory_allocation() {
    local file="$1"
    local file_size=$(du -m "$file" | cut -f1)  # Get size in MB

    if ((file_size < 100)); then
        echo "2G"
    elif ((file_size < 500)); then
        echo "4G"
    elif ((file_size < 1000)); then
        echo "8G"
    else
        echo "16G"
    fi
}

###############################
#  JOB SUBMISSION LOOP
###############################

INDEX=0
while [[ $INDEX -lt $TOTAL_FILES ]]; do
    SAMPLE="${FILES_TO_PROCESS[$INDEX]}"

    # Determine memory requirement based on input file size
    MEM_REQUIRED=$(get_memory_allocation "$INPUT_DIR/${SAMPLE}_1.fastq.gz")

    ###############################
    #  CHECK CURRENT JOB COUNT
    ###############################
    while true; do
        RUNNING_JOBS=$(squeue -u "$USER" -o "%A %j" | grep -c "fastq_")

        if [[ $RUNNING_JOBS -lt $MAX_JOBS ]]; then
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
#SBATCH --time=02:00:00
#SBATCH --mem=${MEM_REQUIRED}
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

# Directories and tools
INPUT_DIR="$INPUT_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
CHECKPOINT_FILE="$CHECKPOINT_FILE"
ADAPTER_FILE="$ADAPTER_FILE"  # Built-in adapters
THREADS=8
SAMPLE_NAME="$SAMPLE"

# Potential paired-end files
R1="\$INPUT_DIR/\${SAMPLE_NAME}_1.fastq.gz"
R2="\$INPUT_DIR/\${SAMPLE_NAME}_2.fastq.gz"

if [[ -f "\$R1" && -f "\$R2" ]]; then
    echo "Processing PAIRED-END sample: \$SAMPLE_NAME"

    # BBDuk for paired-end (adapter trim, mask Ns, minlen)
    bbduk.sh \\
        in1="\$R1" \\
        in2="\$R2" \\
        out1="\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_1.fastq.gz" \\
        out2="\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed_2.fastq.gz" \\
        ref="\$ADAPTER_FILE" \\
        ktrim=r \\          # Trim adapters from the right
        k=23 \\             # Kmer length for adapter matching
        mink=11 \\          # Allow shorter kmers at read ends
        hdist=1 \\          # Hamming distance for matching
        tpe \\              # Trim paired-end reads to same length
        tbo \\              # Trim by overlap detection
        maq=20 \\           # Mask bases with quality <20 as N
        minlen=50 \\        # Discard reads shorter than 50bp
        threads=\$THREADS

elif [[ -f "\$INPUT_DIR/\${SAMPLE_NAME}.fastq.gz" ]]; then
    echo "Processing SINGLE-END sample: \$SAMPLE_NAME"

    # BBDuk for single-end (adapter trim, mask Ns, minlen)
    bbduk.sh \\
        in="\$INPUT_DIR/\${SAMPLE_NAME}.fastq.gz" \\
        out="\$OUTPUT_DIR/\${SAMPLE_NAME}_trimmed.fastq.gz" \\
        ref="\$ADAPTER_FILE" \\
        ktrim=r \\
        k=23 \\
        mink=11 \\
        hdist=1 \\
        maq=20 \\
        minlen=50 \\
        threads=\$THREADS

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
    echo "Submitted job for sample: $SAMPLE with memory $MEM_REQUIRED"

    # Move to the next sample
    ((INDEX++))
done

echo "Done submitting jobs."
