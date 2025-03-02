#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=filteredreads/hmp_parallel_filter.log
#SBATCH --error=filteredreads/hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low

# Define base and output directories
BASE_DIR="hmp_16s_trimmed"
OUTPUT_DIR="filteredreads"

# Create output and logs directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Find all relevant files (.fasta.bz2, .fsa.bz2, .fastq.bz2) and read them into an array
mapfile -t files < <(find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" \))
TOTAL_FILES=${#files[@]}
echo "Total files to process: $TOTAL_FILES"

# Function to determine memory allocation based on file size (in MB)
get_memory_allocation() {
    local file="$1"
    local file_size
    file_size=$(du -m "$file" | cut -f1)
    if (( file_size < 100 )); then
        echo "2G"
    elif (( file_size < 500 )); then
        echo "4G"
    elif (( file_size < 1000 )); then
        echo "8G"
    else
        echo "16G"
    fi
}

# Function to construct the filtering command based on file type.
# It outputs the filtered file to the OUTPUT_DIR while keeping the original filename.
get_filter_command() {
    local file="$1"
    local base_name
    base_name=$(basename "$file")
    local output_file="$OUTPUT_DIR/$base_name"

    if [[ "$file" == *.fasta.bz2 || "$file" == *.fsa.bz2 ]]; then
        # For FASTA files: skip the first empty record and filter records where the first sequence line is at least 30 bp.
        echo "bzcat \"$file\" | awk 'BEGIN {RS=\">\"; ORS=\"\"} NR>1 {if(length(\$2) >= 30) print \">\"\$0}' | bzip2 > \"$output_file\""
    elif [[ "$file" == *.fastq.bz2 ]]; then
        # For FASTQ files: process four lines at a time; only print the full record if the sequence line (second line) is at least 30 bp.
        echo "bzcat \"$file\" | awk '{
            header=\$0;
            getline seq;
            getline plus;
            getline qual;
            if(length(seq) >= 30) {
                print header;
                print seq;
                print plus;
                print qual;
            }
        }' | bzip2 > \"$output_file\""
    fi
}

# Function to submit a job with dynamic memory allocation for a single file
submit_job() {
    local file="$1"
    local mem_required
    mem_required=$(get_memory_allocation "$file")
    local filter_command
    filter_command=$(get_filter_command "$file")

    sbatch --job-name=filter_job \
           --output=logs/%j.out \
           --error=logs/%j.err \
           --cpus-per-task=1 \
           --mem="$mem_required" \
           --time=4:00:00 \
           --wrap="$filter_command"
}

# Maximum number of concurrent jobs allowed
MAX_JOBS=50

# Submit jobs for each file while ensuring we never exceed MAX_JOBS running concurrently
for file in "${files[@]}"; do
    # Wait if the number of running filter jobs has reached MAX_JOBS
    while (( $(squeue --noheader -n filter_job | wc -l) >= MAX_JOBS )); do
        sleep 10
    done
    submit_job "$file"
done

echo "All jobs submitted. Waiting for completion..."

# Wait until all submitted filter jobs are finished
while squeue --noheader -n filter_job | grep -q '[0-9]'; do
    sleep 60
    echo "Waiting for remaining jobs to finish..."
done

echo "All jobs are complete!"
