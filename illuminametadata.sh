#!/bin/bash
#SBATCH --job-name=fastq_processing
#SBATCH --output=fastq_processing_%j.out
#SBATCH --error=fastq_processing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  # Adjust based on your cluster
#SBATCH --mem=16G           # Adjust as needed
#SBATCH --time=48:00:00     # Adjust time limit (hh:mm:ss)

# Load any necessary modules (if required, e.g., for parallel)
module load parallel

# Directory containing FASTQ files
INPUT_DIR="./fastq_files"
OUTPUT_FILE="output.tsv"
LOCKFILE="/tmp/fastq_processing.lock"

# Ensure output file exists with header row
echo -e "Filename\tAccession\tInstrument\tRun\tFlowcell_ID\tLane\tTile\tX\tY\tRead\tIsFiltered\tControl\tIndex" > "$OUTPUT_FILE"

# Function to process a single FASTQ file
process_fastq() {
    file="$1"
    accession=$(basename "$file" .fastq.gz)  # Extract accession

    # Extract header lines, parse, and append to the output file
    zcat "$file" | sed -n '1~4p' | awk -v filename="$file" -v accession="$accession" '
    BEGIN { OFS="\t" }
    {
        match($0, /^@([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^ ]+) ([^:]+):([^:]+):([^:]+):([^ ]+)/, a);
        if (RLENGTH > 0) {
            print filename, accession, a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11];
        }
    }' | flock "$LOCKFILE" tee -a "$OUTPUT_FILE" >/dev/null
}

export -f process_fastq
export OUTPUT_FILE LOCKFILE

# Get number of available CPUs from SLURM
NUM_CPUS=$SLURM_CPUS_PER_TASK
if [ -z "$NUM_CPUS" ]; then
    NUM_CPUS=4  # Default to 4 if not running in SLURM
fi

echo "Using $NUM_CPUS CPUs for parallel processing"

# Run in parallel across all FASTQ files in the directory
find "$INPUT_DIR" -type f -name "*.fastq.gz" | parallel --will-cite -j "$NUM_CPUS" process_fastq {}

echo "Processing complete. Output saved to $OUTPUT_FILE."
