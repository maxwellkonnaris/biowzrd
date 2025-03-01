#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=hmp_parallel_filter.log
#SBATCH --error=hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=low
#SBATCH --array=1-50%50

# Set the base directory containing the hmp_16s files
BASE_DIR="hmp_16s_trimmed"

# Find all .fa.bz2 files and store them in a temporary list
FILE_LIST=$(mktemp)
find "$BASE_DIR" -type f -name "*.fa.bz2" > "$FILE_LIST"
TOTAL_FILES=$(wc -l < "$FILE_LIST")

# Function to determine memory based on file size
get_memory_allocation() {
    local file_size
    file_size=$(du -m "$1" | cut -f1)  # Get file size in MB

    if ((file_size < 100)); then
        echo "2G"
    elif ((file_size < 500)); then
        echo "4G"
    else
        echo "8G"
    fi
}

# Function to submit jobs with dynamic memory allocation
submit_job() {
    local file="$1"
    local mem_required
    mem_required=$(get_memory_allocation "$file")

    sbatch --job-name=filter_job \
           --output=logs/%j.out \
           --error=logs/%j.err \
           --cpus-per-task=1 \
           --mem="$mem_required" \
           --time=4:00:00 \
           --wrap="
    bzcat \"$file\" | awk 'BEGIN {RS=\">\"; ORS=\"\"} length(\$2) >= 30 {print \">\"\$0}' | bzip2 > \"${file}.tmp.bz2\";
    if [[ -s \"${file}.tmp.bz2\" ]]; then mv \"${file}.tmp.bz2\" \"$file\"; else rm -f \"${file}.tmp.bz2\"; fi
    "
}

# Create a logs directory if it doesn't exist
mkdir -p logs

# Submit initial batch of jobs (up to 50 at a time)
for i in $(seq 1 50); do
    read -r file <&3 || break
    submit_job "$file"
done 3<"$FILE_LIST"

# Wait for jobs to finish and submit new ones as slots open up
while squeue --name=filter_job --format="%i" | grep -q '[0-9]'; do
    sleep 10  # Wait a bit before checking again

    # Check running job count
    running_jobs=$(squeue --name=filter_job --noheader | wc -l)

    # If there are available slots, submit more jobs
    while [[ "$running_jobs" -lt 50 ]]; do
        read -r file <&3 || break
        submit_job "$file"
        running_jobs=$((running_jobs + 1))
    done 3<"$FILE_LIST"
done

echo "All jobs submitted and processed!"

# Cleanup
rm -f "$FILE_LIST"

