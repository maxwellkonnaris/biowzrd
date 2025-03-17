#!/bin/bash
#SBATCH --job-name=process16s
#SBATCH --output=process16s_%j.out
#SBATCH --error=process16s_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

INPUT_FASTQ="$1"
SAMPLE_TYPE="$2"

# Validate input
if [[ -z "$INPUT_FASTQ" || -z "$SAMPLE_TYPE" ]]; then
  echo "Usage: $0 <fastq.gz> <sample_type>"
  echo "  <sample_type>: '16S' for DADA2, 'meta' for MetaPhlAn4 + mOTUv2.5"
  exit 1
fi

# Set resources dynamically
FILE_SIZE_GB=$(du -BG "$INPUT_FASTQ" | cut -f1 | sed 's/G//')

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

# Create the job script
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

EOF

# Add workflow-specific commands
if [[ "$SAMPLE_TYPE" == "16S" ]]; then
  cat <<EOF >> "$JOB_SCRIPT"
# ----------------------------
# 16S workflow (DADA2)
# ----------------------------
conda run -n dada2 fastp -i "$INPUT_FASTQ" -o "qc_${INPUT_FASTQ}"
conda run -n dada2 Rscript run_dada2.R "qc_${INPUT_FASTQ}"
EOF
elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
  cat <<EOF >> "$JOB_SCRIPT"
# ----------------------------
# Metagenomics workflow
# ----------------------------
# QC + MetaPhlAn4
conda run -n metaphlan fastp -i "$INPUT_FASTQ" -o "qc_${INPUT_FASTQ}"
conda run -n metaphlan metaphlan "qc_${INPUT_FASTQ}" --input_type fastq --output "metaphlan_${INPUT_FASTQ}.txt"

# mOTUs (reuse QC'd file)
conda run -n motus motus profile -s "qc_${INPUT_FASTQ}" -o "motus_${INPUT_FASTQ}.txt"
EOF
else
  echo "Invalid sample type: $SAMPLE_TYPE"
  exit 1
fi

# Submit the job
sbatch "$JOB_SCRIPT"
echo "Job submitted: $JOB_SCRIPT"
