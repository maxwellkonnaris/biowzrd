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

# Extract sample name and find paired file
BASENAME=$(basename "$INPUT_FASTQ" | sed -E 's/_[12]\.fastq\.gz//')
DIRNAME=$(dirname "$INPUT_FASTQ")

if [[ "$INPUT_FASTQ" =~ "_1.fastq.gz" ]]; then
  PAIRED_FASTQ="${DIRNAME}/${BASENAME}_2.fastq.gz"
elif [[ "$INPUT_FASTQ" =~ "_2.fastq.gz" ]]; then
  PAIRED_FASTQ="${DIRNAME}/${BASENAME}_1.fastq.gz"
else
  PAIRED_FASTQ=""
fi

# Check if the paired file exists and prevent duplicate processing
LOG_FILE="processed_files.log"

if [[ -f "$PAIRED_FASTQ" ]]; then
  if grep -q "$PAIRED_FASTQ" "$LOG_FILE" || grep -q "$INPUT_FASTQ" "$LOG_FILE"; then
    echo "Skipping $INPUT_FASTQ as its pair ($PAIRED_FASTQ) has already been processed."
    exit 0
  fi
  echo "$INPUT_FASTQ" >> "$LOG_FILE"
  echo "$PAIRED_FASTQ" >> "$LOG_FILE"
else
  if grep -q "$INPUT_FASTQ" "$LOG_FILE"; then
    echo "Skipping $INPUT_FASTQ as it has already been processed."
    exit 0
  fi
  echo "$INPUT_FASTQ" >> "$LOG_FILE"
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
JOB_SCRIPT="job_${SAMPLE_TYPE}_${BASENAME}.sh"

# Export global variables
export INPUT_FASTQ="$INPUT_FASTQ"
export SAMPLE_TYPE="$SAMPLE_TYPE"
export PAIRED_FASTQ="$PAIRED_FASTQ"
export BASENAME="$BASENAME"
export JOB_SCRIPT="$JOB_SCRIPT"
export FILE_SIZE_GB="$FILE_SIZE_GB"
export CPUS="$CPUS"
export MEM="$MEM"
export TIME="$TIME"

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

# Export global variables
export INPUT_FASTQ="$INPUT_FASTQ"
export SAMPLE_TYPE="$SAMPLE_TYPE"
export PAIRED_FASTQ="$PAIRED_FASTQ"
export BASENAME="$BASENAME"
export JOB_SCRIPT="$JOB_SCRIPT"
export FILE_SIZE_GB="$FILE_SIZE_GB"
export CPUS="$CPUS"
export MEM="$MEM"
export TIME="$TIME"

EOF

# Add workflow-specific commands
if [[ "$SAMPLE_TYPE" == "16S" ]]; then
  cat <<EOF >> "$JOB_SCRIPT"
# ----------------------------
# 16S workflow (DADA2)
# ----------------------------
if [[ -f "$PAIRED_FASTQ" ]]; then
  echo "Processing paired-end reads: $INPUT_FASTQ and $PAIRED_FASTQ"
  conda run -n dada2 fastp -i "$INPUT_FASTQ" -I "$PAIRED_FASTQ" -o "qc_${BASENAME}_1.fastq.gz" -O "qc_${BASENAME}_2.fastq.gz"
  conda run -n dada2 Rscript run_dada2.R "qc_${BASENAME}_1.fastq.gz" "qc_${BASENAME}_2.fastq.gz"
else
  echo "Processing single-end read: $INPUT_FASTQ"
  conda run -n dada2 fastp -i "$INPUT_FASTQ" -o "qc_${INPUT_FASTQ}"
  conda run -n dada2 Rscript run_dada2.R "qc_${INPUT_FASTQ}"
fi
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

# Make job script executable
chmod +x "$JOB_SCRIPT"

# Submit the job
sbatch "$JOB_SCRIPT"
echo "Job submitted: $JOB_SCRIPT"
