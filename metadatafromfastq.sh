#!/bin/bash
#SBATCH --job-name=fastq_processing
#SBATCH --output=fastq_processing_%j.out
#SBATCH --error=fastq_processing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=48:00:00

module load parallel

# Default values
INPUT_DIR="./fastq_files"
OUTPUT_FILE="output.tsv"

# Parse command-line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_DIR="$2"; shift ;;
        -o|--output) OUTPUT_FILE="$2"; shift ;;
        -h|--help)
            echo "Usage: $0 [-i INPUT_DIR] [-o OUTPUT_FILE]"
            exit 0 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

LOCKFILE="/tmp/fastq_processing.lock"

# Ensure output file exists with header row
echo -e "Filename\tAccession\tRead\tFormat\tInstrument\tRun\tFlowcell_ID\tLane\tTile\tX\tY\tReadLength\tIsFiltered\tControl\tIndex" > "$OUTPUT_FILE"

# Function to process a single FASTQ file
process_fastq() {
    file="$1"
    accession=$(basename "$file" | grep -oE '^(SRR|ERR|DRR)[0-9]+')  # Extract correct SRA accession

    # Extract header lines, parse, and append to the output file
    zcat "$file" | sed -n '1~4p' | awk -v filename="$file" -v accession="$accession" '
    BEGIN { OFS="\t" }
    {
        # 1. SRA/ENA/DDBJ Format (e.g., @SRR27930322.1 M03047:101:000000000-KMJW5:1:1101:16949:1000 length=246)
        if ($1 ~ /^@(SRR|ERR|DRR)/) {
            match($0, /^@([^ ]+)\.([0-9]+) ([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^ ]+) length=([0-9]+)/, a);
            if (RLENGTH > 0) {
                print filename, accession, a[2], "SRA/ENA/DDBJ", a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], "-", "-", "-";
            }
        
        # 2. Illumina CASAVA 1.8+ Format
        } else if ($1 ~ /^@[A-Za-z0-9]+:[0-9]+:/) {
            match($0, /^@([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^ ]+) ([^:]+):([^:]+):([^:]+):([^ ]+)/, a);
            if (RLENGTH > 0) {
                print filename, accession, "-", "Illumina-CASAVA1.8+", a[1], a[2], a[3], a[4], a[5], a[6], a[7], "-", a[8], a[9], a[10];
            }
        
        # 3. Old Illumina Format
        } else if ($1 ~ /^@[A-Za-z0-9-]+:[0-9]+:/) {
            match($0, /^@([^:]+):([^:]+):([^:]+):([^:]+):([^ ]+)#([^ ]+)\/([12])/, a);
            if (RLENGTH > 0) {
                print filename, accession, a[7], "Illumina-Old", a[1], "-", "-", a[2], a[3], a[4], a[5], "-", "-", "-", a[6];
            }

        # 4. Restored SRA Format
        } else if ($1 ~ /^@[0-9]+_[A-Za-z0-9-]+_s_[0-9]+:/) {
            match($0, /^@([^_]+)_([^_]+)_s_([^:]+):([^:]+):([^:]+):([^:]+):([^ ]+)/, a);
            if (RLENGTH > 0) {
                print filename, accession, "-", "SRA-Restored", "-", "-", "-", a[3], a[4], a[5], a[6], "-", "-", "-", "-";
            }
        }
    }' | flock "$LOCKFILE" tee -a "$OUTPUT_FILE" >/dev/null
}

export -f process_fastq
export OUTPUT_FILE LOCKFILE

# Get number of available CPUs from SLURM
NUM_CPUS=$SLURM_CPUS_PER_TASK
if [ -z "$NUM_CPUS" ]; then
    NUM_CPUS=4
fi

echo "Using $NUM_CPUS CPUs for parallel processing"

# Run in parallel across all FASTQ files in the directory
find "$INPUT_DIR" -type f -name "*.fastq.gz" | parallel --will-cite -j "$NUM_CPUS" process_fastq {}

echo "Processing complete. Output saved to $OUTPUT_FILE."
