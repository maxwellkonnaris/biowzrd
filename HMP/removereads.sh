#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=logs/hmp_parallel_filter.log
#SBATCH --error=logs/hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Directories
BASE_DIR="hmp_16s_trimmed"  # Change to your input directory
OUTPUT_DIR="filteredreads"   # Change to your desired output directory
TMP_DIR="temp_processing"    # Temporary directory to ensure no overwrites

mkdir -p "$OUTPUT_DIR" "$TMP_DIR" logs

# Find all relevant files
mapfile -t files < <(find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" -o -name "*.fq.bz2" \))
echo "$(date): Total files to process: ${#files[@]}" | tee -a logs/hmp_parallel_filter.out

# Max concurrent jobs
MAX_JOBS=50

# Function to submit filtering job
submit_job() {
    local file="$1"
    local base_name
    base_name=$(basename "$file")
    local temp_file="$TMP_DIR/$base_name"
    local output_file="$OUTPUT_DIR/$base_name"

    echo "$(date): Processing $file" >> logs/hmp_parallel_filter.out

    # Check if the file exists and is readable
    if [[ ! -f "$file" ]]; then
        echo "$(date): ERROR - File not found: $file" >> logs/hmp_parallel_filter.err
        return
    elif [[ ! -r "$file" ]]; then
        echo "$(date): ERROR - Cannot read file: $file" >> logs/hmp_parallel_filter.err
        return
    fi

    # Copy file to temp directory before processing (preserves original)
    cp "$file" "$temp_file"

    # Verify `bzcat` produces output
    if ! bzcat "$temp_file" | head -n 1 &>/dev/null; then
        echo "$(date): ERROR - bzcat failed or file is empty: $file" >> logs/hmp_parallel_filter.err
        return
    fi

    # Determine filtering logic based on file type
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

    # Submit Slurm job and append logs instead of overwriting
    sbatch \
      --job-name=filter_job \
      --cpus-per-task=1 \
      --mem=2G \
      --time=4:00:00 \
      --wrap="$filter_cmd" >> logs/hmp_parallel_filter.out 2>> logs/hmp_parallel_filter.err

    # Ensure the filtered file is not unexpectedly empty
    sleep 2  # Allow time for file creation
    if [[ -f "$output_file" && $(stat -c%s "$output_file") -lt 100 ]]; then
        echo "$(date): WARNING - Output file too small: $output_file" >> logs/hmp_parallel_filter.err
    fi
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
