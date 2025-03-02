#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=logs/hmp_parallel_filter.log
#SBATCH --error=logs/hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Directories
BASE_DIR="hmp_16s_trimmed"
OUTPUT_DIR="filteredreads"
TMP_DIR="temp_processing"

mkdir -p "$OUTPUT_DIR" "$TMP_DIR" logs

# Find all relevant files
mapfile -t files < <(find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" -o -name "*.fq.bz2" \) -print | sort)

echo "$(date): Total files to process: ${#files[@]}" | tee -a logs/hmp_parallel_filter.out

# Max concurrent jobs
MAX_JOBS=50

# Function to submit filtering job
submit_job() {
    local file="$1"
    local base_name
    base_name=$(basename "$file")
    local temp_file="$TMP_DIR/${base_name}_$$.bz2"  # Unique temp file
    local output_file="$OUTPUT_DIR/$base_name"

    echo "$(date): Checking file existence: $file" >> logs/hmp_parallel_filter.out

    # Verify file exists before processing
    if [[ ! -f "$file" ]]; then
        echo "$(date): ERROR - File not found: $file" >> logs/hmp_parallel_filter.err
        return
    fi

    # Copy file to a unique temp file to prevent overwrite issues
    cp "$file" "$temp_file"

    # Verify `bzcat` can read the file
    if ! bzcat "$temp_file" | head -n 1 &>/dev/null; then
        echo "$(date): ERROR - bzcat failed or empty file: $file" >> logs/hmp_parallel_filter.err
        return
    fi

    echo "$(date): Processing $file" >> logs/hmp_parallel_filter.out

    # Check file type and determine filtering command
    if [[ "$file" == *.fasta.bz2 || "$file" == *.fsa.bz2 ]]; then
        filter_cmd="bzcat \"$temp_file\" | awk '{
            if (\$0 ~ /^>/) {header = \$0; next}
            seq = \$0;
            if (length(seq) >= 30) {
                print header;
                print seq;
            }
        }' | bzip2 > \"$output_file\""

    elif [[ "$file" == *.fastq.bz2 || "$file" == *.fq.bz2 ]]; then
        filter_cmd="bzcat \"$temp_file\" | awk '{
            if(NR%4 == 1) { h=\$0 }
            else if(NR%4 == 2) { s=\$0 }
            else if(NR%4 == 3) { p=\$0 }
            else if(NR%4 == 0) { 
                q=\$0;
                if(length(s) >= 30) {
                    print h;
                    print s;
                    print p;
                    print q;
                }
            }
        }' | bzip2 > \"$output_file\""

    else
        echo "$(date): Skipping unsupported file type: $file" >> logs/hmp_parallel_filter.out
        return
    fi

    # Submit Slurm job and append logs
    sbatch \
      --job-name=filter_job \
      --cpus-per-task=1 \
      --mem=2G \
      --time=4:00:00 \
      --wrap="$filter_cmd" >> logs/hmp_parallel_filter.out 2>> logs/hmp_parallel_filter.err

    # Wait for the output file to be fully written
    sleep 5  # Give the system some time to complete writes

    # Ensure the file exists before proceeding
    WAIT_TIME=0
    MAX_WAIT=60  # Max wait time (seconds)

    while [[ ! -f "$output_file" || $(stat -c%s "$output_file") -lt 100 ]]; do
        if [[ $WAIT_TIME -ge $MAX_WAIT ]]; then
            echo "$(date): ERROR - Output file too small even after waiting: $output_file" >> logs/hmp_parallel_filter.err
            break
        fi
        sleep 2
        WAIT_TIME=$((WAIT_TIME + 2))
    done

    echo "$(date): Verified write completed for $output_file" >> logs/hmp_parallel_filter.out

    # Cleanup: Remove temp file after job submission
    rm -f "$temp_file"
}

# Process files while respecting MAX_JOBS
for file in "${files[@]}"; do
    while (( $(squeue --noheader --format=%j | grep -c '^filter_job$') >= MAX_JOBS )); do
        sleep 10
    done
    submit_job "$file"
done

echo "$(date): All filter jobs submitted. Waiting for completion..." | tee -a logs/hmp_parallel_filter.out

# Wait until all jobs finish
while squeue --noheader --format=%j | grep -q '^filter_job$'; do
    sleep 60
    echo "$(date): Waiting for remaining filter jobs to finish..." | tee -a logs/hmp_parallel_filter.out
done

echo "$(date): All jobs are complete!" | tee -a logs/hmp_parallel_filter.out
