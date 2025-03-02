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

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Find all relevant files
mapfile -t files < <(find "$BASE_DIR" -type f \( -name "*.fasta.bz2" -o -name "*.fsa.bz2" -o -name "*.fastq.bz2" -o -name "*.fq.bz2" \))
echo "Total files to process: ${#files[@]}"

# Max concurrent jobs
MAX_JOBS=50

# Function to submit filtering job
submit_job() {
    local file="$1"
    local base_name
    base_name=$(basename "$file")
    local output_file="$OUTPUT_DIR/$base_name"

    # Check file type and build the filtering command
    if [[ "$file" == *.fasta.bz2 || "$file" == *.fsa.bz2 ]]; then
        # Process FASTA files (2-line records)
        filter_cmd="bzcat \"$file\" | awk '{
            if (\$0 ~ /^>/) {header = \$0; next} 
            seq = \$0;
            if (length(seq) >= 30) {
                print header;
                print seq;
            }
        }' | bzip2 > \"$output_file\""

    elif [[ "$file" == *.fastq.bz2 || "$file" == *.fq.bz2 ]]; then
        # Process FASTQ files (4-line records)
        filter_cmd="bzcat \"$file\" | awk '{
            if(NR%4 == 1) { h=\$0 }         # Header line
            else if(NR%4 == 2) { s=\$0 }    # Sequence line
            else if(NR%4 == 3) { p=\$0 }    # Plus sign (+)
            else if(NR%4 == 0) { 
                q=\$0;                      # Quality line
                if(length(s) >= 30) {       # Keep only if sequence length >= 30
                    print h;
                    print s;
                    print p;
                    print q;
                }
            }
        }' | bzip2 > \"$output_file\""

    else
        echo "Skipping unsupported file type: $file"
        return
    fi

    # Submit as a Slurm job
    sbatch \
      --job-name=filter_job \
      --output=logs/%j.out \
      --error=logs/%j.err \
      --cpus-per-task=1 \
      --mem=2G \
      --time=4:00:00 \
      --wrap="$filter_cmd"
}

# Process files while respecting MAX_JOBS
for file in "${files[@]}"; do
    while (( $(squeue --noheader --format=%j | grep -c '^filter_job$') >= MAX_JOBS )); do
        sleep 10
    done
    submit_job "$file"
done

echo "All filter jobs submitted. Waiting for completion..."

# Wait until all jobs finish
while squeue --noheader --format=%j | grep -q '^filter_job$'; do
    sleep 60
    echo "Waiting for remaining filter jobs to finish..."
done

echo "All jobs are complete!"
