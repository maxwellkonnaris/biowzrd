#!/bin/bash
#SBATCH --job-name=hmp_parallel_filter
#SBATCH --output=filteredreads/hmp_parallel_filter.log
#SBATCH --error=filteredreads/hmp_parallel_filter.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Adjust paths as needed
BASE_DIR="hmp_16s_trimmed"
OUTPUT_DIR="filteredreads"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Find all relevant .fasta.bz2, .fsa.bz2, .fastq.bz2 files
mapfile -t files < <(find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" \))
echo "Total files to process: ${#files[@]}"

# MAX_JOBS: how many jobs to run concurrently
MAX_JOBS=50

# Helper function to submit Slurm job for a single file
submit_job() {
    local file="$1"

    local base_name
    base_name=$(basename "$file")
    local output_file="$OUTPUT_DIR/$base_name"

    # Build the filtering command depending on FASTA vs FASTQ
    if [[ "$file" == *.fasta.bz2 || "$file" == *.fsa.bz2 ]]; then
        # Each record is exactly 2 lines: header + 1 read line
        # Keep if sequence length ≥ 30
        filter_cmd="bzcat \"$file\" | paste - - | awk -F\"\t\" '{
            header=\$1; seq=\$2;
            if(length(seq) >= 30) {
                print header;
                print seq;
            }
        }' | bzip2 > \"$output_file\""

    elif [[ "$file" == *.fastq.bz2 ]]; then
        # Standard FASTQ: 4 lines per record
        # Keep if sequence line (line #2) ≥ 30
        filter_cmd="bzcat \"$file\" | awk '{
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
        echo "Skipping file with unknown extension: $file"
        return
    fi

    # Submit as a separate Slurm job
    sbatch \
      --job-name=filter_job \
      --output=logs/%j.out \
      --error=logs/%j.err \
      --cpus-per-task=1 \
      --mem=2G \
      --time=4:00:00 \
      --wrap="$filter_cmd"
}

# Loop over files, respecting the MAX_JOBS concurrency limit
for file in "${files[@]}"; do
    # If there are already MAX_JOBS filter_job’s running, wait
    while (( $(squeue --noheader -n filter_job | wc -l) >= MAX_JOBS )); do
        sleep 10
    done
    submit_job "$file"
done

echo "All filter jobs submitted. Waiting for completion..."

# Optionally wait until all filter jobs are done
while squeue --noheader -n filter_job | grep -q '[0-9]'; do
    sleep 60
    echo "Waiting for remaining filter jobs to finish..."
done

echo "All jobs are complete!"
