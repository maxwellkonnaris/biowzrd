#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=hmp_parallel_filter.log
#SBATCH --error=hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low

# Set the base directory containing all subdirectories
BASE_DIR="hmp_16s_trimmed"

# Find all relevant files (.fasta.bz2, .fsa.bz2, .fastq.bz2) in all subdirectories
FILE_LIST=$(mktemp)
find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" \) > "$FILE_LIST"
TOTAL_FILES=$(wc -l < "$FILE_LIST")

echo "Total files to process: $TOTAL_FILES"

# Function to determine memory allocation based on file size
get_memory_allocation() {
    local file_size
    file_size=$(du -m "$1" | cut -f1)  # Get file size in MB

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

# Function to determine file type and filtering method
get_filter_command() {
    local file="$1"

    # Detect if file is FASTA (fasta, fsa) or FASTQ
    if [[ "$file" == *.fasta.bz2 || "$file" == *.fsa.bz2 ]]; then
        echo "bzcat \"$file\" | awk 'BEGIN {RS=\">\"; ORS=\"\"} length(\$2) >= 30 {print \">\"\$0}' | bzip2 > \"${file}.tmp.bz2\""
    elif [[ "$file" == *.fastq.bz2 ]]; then
        echo "bzcat \"$file\" | awk 'NR%4==1 || NR%4==2 && length(\$0) >= 30 || NR%4==3 || NR%4==0' | bzip2 > \"${file}.tmp.bz2\""
    fi
}

# Function to submit a job with dynamic memory allocation
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
           --wrap="
    $filter_command;
    if [[ -s \"${file}.tmp.bz2\" ]]; then mv \"${file}.tmp.bz2\" \"$file\"; else rm -f \"${file}.tmp.bz2\"; fi
    "
}

# Create a logs directory if it doesn't exist
mkdir -p logs

# Counter for submitted jobs
SUBMITTED_JOBS=0

# Submit the first batch of 50 jobs
while [[ "$SUBMITTED_JOBS" -lt 50 && -s "$FILE_LIST" ]]; do
    read -r file <&3 || break
    submit_job "$file"
    ((SUBMITTED_JOBS++))
done 3<"$FILE_LIST"

# Process remaining files dynamically
while [[ -s "$FILE_LIST" ]]; do
    sleep 10  # Wait a bit before checking for open slots

    # Check running job count
    running_jobs=$(squeue --name=filter_job --noheader | wc -l)

    # Submit new jobs as slots free up
    while [[ "$running_jobs" -lt 50 && -s "$FILE_LIST" ]]; do
        read -r file <&3 || break
        submit_job "$file"
        ((running_jobs++))
    done 3<"$FILE_LIST"
done

echo "All jobs submitted. Waiting for completion..."

# Wait for all jobs to finish before exiting
while squeue --name=filter_job --format="%i" | grep -q '[0-9]'; do
    sleep 60
    echo "Waiting for remaining jobs to finish..."
done

echo "All jobs are complete!"
rm -f "$FILE_LIST"

