#!/bin/bash
#SBATCH --job-name=process16s
#SBATCH --output=process16s_%j.out   # Standard output log
#SBATCH --error=process16s_%j.err    # Error log
#SBATCH --time=48:00:00              # Maximum time (adjusted dynamically)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=1            # Default CPU cores (adjusted dynamically)
#SBATCH --mem=4G                     # Default memory (adjusted dynamically)
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu  # Email notifications
#SBATCH --mail-type=END,FAIL         # Send mail on end or failure

# Usage: ./submit_jobs.sh <fastq.gz> <sample_type>
# <sample_type>: "16S" for DADA2 or "meta" for metagenomic

INPUT_FASTQ="$1"
SAMPLE_TYPE="$2"  # "16S" or "meta"

if [[ -z "$INPUT_FASTQ" || -z "$SAMPLE_TYPE" ]]; then
    echo "Usage: $0 <fastq.gz> <sample_type>"
    echo "  <sample_type>: '16S' for DADA2, 'meta' for MetaPhlAn4 + mOTUv2.5"
    exit 1
fi

# Determine file size in GB
FILE_SIZE_GB=$(du -BG "$INPUT_FASTQ" | cut -f1 | sed 's/G//')

# Set resource allocation based on file size
if (( FILE_SIZE_GB <= 5 )); then
    CPUS=2
    MEM="8G"
    TIME="12:00:00"
elif (( FILE_SIZE_GB <= 20 )); then
    CPUS=4
    MEM="16G"
    TIME="24:00:00"
else
    CPUS=8
    MEM="32G"
    TIME="48:00:00"
fi

echo "Submitting job for $INPUT_FASTQ with $CPUS CPUs, $MEM memory, and $TIME runtime..."

# Create the SLURM job script dynamically
JOB_SCRIPT="job_${SAMPLE_TYPE}_${SLURM_JOB_ID}.sh"

cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH --job-name=process_${SAMPLE_TYPE}
#SBATCH --output=process_${SAMPLE_TYPE}_%j.out
#SBATCH --error=process_${SAMPLE_TYPE}_%j.err
#SBATCH --time=$TIME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

# Quality control using fastp
fastp -i "$INPUT_FASTQ" -o "qc_${INPUT_FASTQ}"

EOF

# Add specific processing steps
if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    cat <<EOF >> "$JOB_SCRIPT"
# Running DADA2 for 16S data
Rscript run_dada2.R "qc_${INPUT_FASTQ}"
EOF
elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    cat <<EOF >> "$JOB_SCRIPT"
# Running MetaPhlAn4 for taxonomic profiling
metaphlan "qc_${INPUT_FASTQ}" --input_type fastq --output "metaphlan_${INPUT_FASTQ}.txt"

# Running mOTUv2.5 for functional profiling
motus profile -s "qc_${INPUT_FASTQ}" -o "motus_${INPUT_FASTQ}.txt"
EOF
else
    echo "Error: Invalid sample type. Use '16S' or 'meta'."
    exit 1
fi

# Submit the job
sbatch "$JOB_SCRIPT"

echo "Job submitted: $JOB_SCRIPT"
