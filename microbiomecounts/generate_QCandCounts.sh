#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --output=counts_%j.out
#SBATCH --error=counts_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

# Default directory for FASTQ files
DEFAULT_DIR="fastq_data/fastq_biologicaldata/"
LOCK_DIR=$(mktemp -d)
DEBUG_FILE="debug.log"
DEBUG_LOCK="$LOCK_DIR/debug.lock"
COMPLETED_FILE="completed_steps.log"
COMPLETED_LOCK="$LOCK_DIR/completed.lock"
FAILED_FILE="failed.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
DEFAULT_WORKERS=8  # Default number of workers
DEFAULT_MOTUS_TAX_LEVEL="mOTU"  # Default taxonomic level for mOTUs
LOG_LEVEL="INFO"  # Adjustable log level: INFO, DEBUG
TMP_BASE="${SCRATCH:-/scratch}"  # Prefer $SCRATCH, fall back to /scratch
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"
QUALITY_PROFILE_DIR="quality_profiles"
# Default micromamba environment paths or names
DEFAULT_DADA2_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/dada2"
DEFAULT_MOTUS_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/motus"
DEFAULT_MPA_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa"

# Cleanup function for lock directory
cleanup() {
  rm -rf "$LOCK_DIR"
}
trap cleanup EXIT

# Parse command-line arguments
while getopts ":i:d:w:k:qa:m:p:" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    d) FASTQ_DIR="$OPTARG" ;;
    w) NUM_WORKERS="$OPTARG" ;;
    k) MOTUS_TAX_LEVEL="$OPTARG" ;;
    q) QUALITY_CHECK="true" ;;
    a) DADA2_ENV="$OPTARG" ;;
    m) MOTUS_ENV="$OPTARG" ;;
    p) MPA_ENV="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Set defaults
FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_DIR}"
NUM_WORKERS="${NUM_WORKERS:-$DEFAULT_WORKERS}"
MOTUS_TAX_LEVEL="${MOTUS_TAX_LEVEL:-$DEFAULT_MOTUS_TAX_LEVEL}"
DADA2_ENV="${DADA2_ENV:-$DEFAULT_DADA2_ENV}"
MOTUS_ENV="${MOTUS_ENV:-$DEFAULT_MOTUS_ENV}"
MPA_ENV="${MPA_ENV:-$DEFAULT_MPA_ENV}"

# Validate micromamba and environments
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found"
  exit 1
fi

# Function to check if an environment exists
check_env() {
  local env="$1"
  local env_name="$2"
  # If env is a path, check if it exists as a directory
  if [[ -d "$env" ]]; then
    return 0
  # Otherwise, assume it's a name and check micromamba env list
  elif micromamba env list | grep -qE "^[[:space:]]*$env[[:space:]]"; then
    return 0
  else
    echo "Error: micromamba environment $env_name ($env) not found"
    exit 1
  fi
}

# Validate each environment
check_env "$DADA2_ENV" "dada2"
check_env "$MOTUS_ENV" "motus"
check_env "$MPA_ENV" "mpa"

# Validate input file and directories
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv|input.csv|input.txt> [-d <fastq_directory>] [-w <num_workers>] [-k <motus_tax_level>] [-q] [-a <dada2_env>] [-m <motus_env>] [-p <mpa_env>]"
  echo "  <input>: CSV, TSV, or TXT file with columns: Bioproject,RunAccession,SequencingType"
  echo "  <fastq_directory>: Directory containing FASTQ files (default: $DEFAULT_DIR)"
  echo "  <num_workers>: Number of parallel workers (default: $DEFAULT_WORKERS)"
  echo "  <motus_tax_level>: Taxonomic level for mOTUs (e.g., mOTU, phylum; default: $DEFAULT_MOTUS_TAX_LEVEL)"
  echo "  -q: Generate quality profiles for FASTQs before processing"
  echo "  <dada2_env>: DADA2 environment path or name (default: $DEFAULT_DADA2_ENV)"
  echo "  <motus_env>: mOTUs environment path or name (default: $DEFAULT_MOTUS_ENV)"
  echo "  <mpa_env>: MetaPhlAn environment path or name (default: $DEFAULT_MPA_ENV)"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file $INPUT_FILE does not exist."
  exit 1
fi

# Validate temporary directory base
if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
  echo "Warning: $TMP_BASE not available or not writable, falling back to /tmp"
  TMP_BASE="/tmp"
  if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
    echo "Error: No writable temporary directory base available (/scratch or /tmp)."
    exit 1
  fi
fi

# Detect input file format
validate_input_file() {
  local file="$1"
  local first_line
  first_line=$(head -n 1 "$file")
  
  # Check for CSV (comma-separated)
  if [[ "$first_line" == *"Bioproject,RunAccession,SequencingType"* ]]; then
    DELIMITER=','
    EXPECTED_HEADER="Bioproject,RunAccession,SequencingType"
  # Check for TSV/TXT (tab-separated)
  elif [[ "$first_line" == $'Bioproject\tRunAccession\tSequencingType'* ]]; then
    DELIMITER=$'\t'
    EXPECTED_HEADER="Bioproject\tRunAccession\tSequencingType"
  else
    echo "Error: Input file $file has unrecognized format or header. Expected: Bioproject,RunAccession,SequencingType (CSV) or Bioproject<tab>RunAccession<tab>SequencingType (TSV/TXT)"
    exit 1
  fi

  if [[ "$first_line" != "$EXPECTED_HEADER" ]]; then
    echo "Error: Input file $file has incorrect header. Expected: $EXPECTED_HEADER"
    exit 1
  fi
}
validate_input_file "$INPUT_FILE"

# Check for 16S samples and validate RDP database
check_rdp_database() {
  local has_16s=false
  while IFS="$DELIMITER" read -r bioproject accession sample_type; do
    if [[ "$sample_type" == "16S" ]]; then
      has_16s=true
      break
    fi
  done < <(tail -n +2 "$INPUT_FILE")
  
  if [[ "$has_16s" == "true" ]]; then
    if [[ ! -f "$RDP_DATABASE" ]]; then
      echo "Error: RDP database $RDP_DATABASE not found in current working directory, required for 16S samples."
      exit 1
    fi
    log_info "RDP database $RDP_DATABASE found for 16S samples."
  fi
}
check_rdp_database

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

# Calculate threads per worker
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
if [[ $THREADS_PER_WORKER -lt 1 ]]; then
  THREADS_PER_WORKER=1
  NUM_WORKERS=$SLURM_CPUS_PER_TASK
fi
echo "Running with $NUM_WORKERS workers, $THREADS_PER_WORKER threads each (total $SLURM_CPUS_PER_TASK CPUs)"

# Logging function with file locking and log level
log_debug() {
  [[ "$LOG_LEVEL" == "DEBUG" ]] || return
  local msg="$1"
  echo "$(date) [DEBUG] $msg" >&2
  flock -x "$DEBUG_LOCK" -c "echo \"$(date) [DEBUG] $msg\" >> \"$DEBUG_FILE\"" 2>/dev/null || \
    echo "[WARN] Could not log to $DEBUG_FILE: $msg" >&2
}

log_info() {
  local msg="$1"
  echo "$(date) [INFO] $msg" >&2
  flock -x "$DEBUG_LOCK" -c "echo \"$(date) [INFO] $msg\" >> \"$DEBUG_FILE\"" 2>/dev/null || \
    echo "[WARN] Could not log to $DEBUG_FILE: $msg" >&2
}

# Append to file with locking
append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  flock -x "$lock" -c "echo \"$line\" >> \"$file\"" 2>/dev/null || \
    echo "[ERROR] Could not append '$line' to $file" >&2
}

# Validate FASTQ file (size and gzip integrity)
validate_fastq() {
  local file="$1"
  local min_size=100  # Minimum size in bytes
  if [[ ! -f "$file" ]]; then
    log_debug "[validate_fastq] File $file does not exist"
    return 1
  fi
  local size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    log_debug "[validate_fastq] File $file too small: $size bytes"
    return 1
  fi
  if ! gzip -t "$file" 2>/dev/null; then
    log_debug "[validate_fastq] File $file is not a valid gzip"
    return 1
  fi
  return 0
}

# Run command with retry
run_command() {
  local cmd="$1"
  local err_msg="$2"
  local max_attempts=3
  local attempt=1
  local wait_time=2

  while (( attempt <= max_attempts )); do
    if eval "$cmd"; then
      return 0
    else
      log_debug "$err_msg (attempt $attempt/$max_attempts)"
      if (( attempt == max_attempts )); then
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

# Run command with retry and capture output
run_command_with_output() {
  local cmd="$1"
  local err_msg="$2"
  local output_file="$3"
  local max_attempts=3
  local attempt=1
  local wait_time=2

  while (( attempt <= max_attempts )); do
    if eval "$cmd > \"$output_file\" 2>&1"; then
      return 0
    else
      log_debug "$err_msg (attempt $attempt/$max_attempts)"
      if (( attempt == max_attempts )); then
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

# Generate quality profiles for FASTQs
generate_quality_profiles() {
  local max_samples=5  # Limit to avoid overwhelming output
  local count=0
  mkdir -p "$QUALITY_PROFILE_DIR"
  while IFS="$DELIMITER" read -r bioproject accession sample_type; do
    if [[ "$sample_type" == "16S" && $count -lt $max_samples ]]; then
      local fastq1="${FASTQ_DIR}/${accession}_1.fastq.gz"
      local fastq2="${FASTQ_DIR}/${accession}_2.fastq.gz"
      local output_pdf="${QUALITY_PROFILE_DIR}/${accession}_quality.pdf"
      if [[ -f "$fastq1" && ! -f "$output_pdf" ]]; then
        if [[ -f "$fastq2" ]]; then
          run_command "micromamba run -n \"$DADA2_ENV\" Rscript -e 'library(dada2); pdf(\"$output_pdf\"); plotQualityProfile(c(\"$fastq1\", \"$fastq2\")); dev.off()'" \
            "[quality profile] Failed for $accession" || log_debug "Failed to generate quality profile for $accession"
        else
          run_command "micromamba run -n \"$DADA2_ENV\" Rscript -e 'library(dada2); pdf(\"$output_pdf\"); plotQualityProfile(\"$fastq1\"); dev.off()'" \
            "[quality profile] Failed for $accession" || log_debug "Failed to generate quality profile for $accession"
        fi
        log_info "Generated quality profile for $accession in $output_pdf"
        (( count++ ))
      fi
    fi
  done < <(tail -n +2 "$INPUT_FILE")
}

# Convert MetaPhlAn profile to counts
convert_metaphlan_to_counts() {
  local metaphlan_log="$1"
  local metaphlan_profile="$2"
  local output_file="$3"
  if [[ ! -f "$metaphlan_log" || ! -f "$metaphlan_profile" ]]; then
    log_debug "[convert_metaphlan_to_counts] Missing files: $metaphlan_log or $metaphlan_profile"
    return 1
  fi
  local mapped_reads=$(grep "Total number of reads mapped" "$metaphlan_log" | awk '{print $6}' | sed 's/(.*//')
  if [[ -z "$mapped_reads" ]]; then
    log_debug "[convert_metaphlan_to_counts] Could not find mapped reads in $metaphlan_log"
    return 1
  fi
  awk -v mapped="$mapped_reads" '
    BEGIN { FS="\t"; OFS="\t"; print "#clade_name\trelative_abundance\tread_count" }
    /^#/ { next }
    { count = ($2 * mapped / 100); print $1, $2, count }
  ' "$metaphlan_profile" > "$output_file"
  log_debug "[convert_metaphlan_to_counts] Converted $metaphlan_profile to counts in $output_file (total mapped: $mapped_reads)"
  return 0
}

# Merge profiles for a Bioproject
merge_profiles() {
  local bioproject="$1"
  local tool="$2"  # "metaphlan", "motus", or "dada2"
  local output_dir="$bioproject"
  local merged_file="${output_dir}/${bioproject}_${tool}_merged.txt"
  local profile_files

  # Load expected accessions for this Bioproject
  declare -A expected_accessions
  while IFS="$DELIMITER" read -r proj accession sample_type; do
    if [[ "$proj" == "$bioproject" ]]; then
      expected_accessions["$accession"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  # Find available profile files
  if [[ "$tool" == "metaphlan" ]]; then
    profile_files=($(ls "$output_dir"/*_metaphlan4_counts.txt 2>/dev/null))
  elif [[ "$tool" == "motus" ]]; then
    profile_files=($(ls "$output_dir"/*_motus.txt 2>/dev/null))
  elif [[ "$tool" == "dada2" ]]; then
    profile_files=($(ls "$output_dir"/seqtab_*.rds 2>/dev/null))
  fi

  # Check if all expected accessions have completed profiles
  local all_complete=true
  for accession in "${!expected_accessions[@]}"; do
    if [[ -z "${COMPLETED_SAMPLES[$accession]}" ]]; then
      log_debug "[merge_profiles] $bioproject ($tool): Accession $accession not fully processed"
      all_complete=false
      break
    fi
    if [[ "$tool" == "metaphlan" && ! -f "$output_dir/${accession}_metaphlan4_counts.txt" ]]; then
      log_debug "[merge_profiles] $bioproject ($tool): Missing metaphlan profile for $accession"
      all_complete=false
      break
    elif [[ "$tool" == "motus" && ! -f "$output_dir/${accession}_motus.txt" ]]; then
      log_debug "[merge_profiles] $bioproject ($tool): Missing motus profile for $accession"
      all_complete=false
      break
    elif [[ "$tool" == "dada2" && ! -f "$output_dir/seqtab_${accession}.rds" ]]; then
      log_debug "[merge_profiles] $bioproject ($tool): Missing dada2 seqtab for $accession"
      all_complete=false
      break
    fi
  done

  # Proceed with merging only if all expected profiles are complete
  if [[ "$all_complete" == "true" && ${#profile_files[@]} -gt 0 ]]; then
    if [[ "$tool" == "metaphlan" && ${#profile_files[@]} -gt 1 ]]; then
      local input_list="${profile_files[*]}"
      run_command "micromamba run -n \"$MPA_ENV\" merge_metaphlan_tables.py $input_list > \"$merged_file\"" \
        "[metaphlan merge] Failed for $bioproject" || return 1
      log_debug "Merged MetaPhlAn profiles for $bioproject into $merged_file"
    elif [[ "$tool" == "metaphlan" && ${#profile_files[@]} -eq 1 ]]; then
      cp "${profile_files[0]}" "$merged_file"
      log_debug "Copied single metaphlan profile for $bioproject to $merged_file"
    elif [[ "$tool" == "motus" && ${#profile_files[@]} -gt 1 ]]; then
      local input_list=$(printf "%s," "${profile_files[@]}" | sed 's/,$//')
      run_command "micromamba run -n \"$MOTUS_ENV\" motus merge -i \"$input_list\" -o \"$merged_file\"" \
        "[motus merge] Failed for $bioproject" || return 1
      log_debug "Merged mOTUs profiles for $bioproject into $merged_file"
    elif [[ "$tool" == "motus" && ${#profile_files[@]} -eq 1 ]]; then
      cp "${profile_files[0]}" "$merged_file"
      log_debug "Copied single motus profile for $bioproject to $merged_file"
    elif [[ "$tool" == "dada2" ]]; then
      local input_list="${bioproject} ${profile_files[*]}"
      run_command "micromamba run -n \"$DADA2_ENV\" Rscript merge_dada2.R $input_list" \
        "[dada2 merge] Failed for $bioproject" || return 1
      log_debug "Merged and processed DADA2 sequence tables for $bioproject"
    fi
  else
    log_debug "[merge_profiles] $bioproject ($tool): Skipping merge due to incomplete or missing profiles"
    return 1
  fi
  return 0
}

# Load completed samples into an associative array
declare -A COMPLETED_SAMPLES
load_completed() {
  if [[ -f "$COMPLETED_FILE" ]]; then
    while IFS=":" read -r accession status; do
      if [[ "$status" == "COMPLETE" ]]; then
        COMPLETED_SAMPLES["$accession"]=1
      fi
    done < "$COMPLETED_FILE"
  fi
  log_debug "Loaded ${#COMPLETED_SAMPLES[@]} completed samples from $COMPLETED_FILE"
}

# Process a single sample
process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SAMPLE_TYPE="$3"

  # Check if already completed
  if [[ -n "${COMPLETED_SAMPLES[$RUN_ACCESSION]}" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already fully processed"
    return
  fi

  # Define file paths
  local INPUT_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_1.fastq.gz"
  local PAIRED_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_2.fastq.gz"
  local OUTPUT_DIR="${BIOPROJECT}"
  local LOG_FILE="${OUTPUT_DIR}/processed_files.log"
  local TMP_DIR=$(mktemp -d -t "process_${SLURM_JOB_ID}_${RUN_ACCESSION}_XXXXXX" -p "$TMP_BASE")
  local QC1="${TMP_DIR}/qc_${RUN_ACCESSION}_1.fastq.gz"
  local QC2="${TMP_DIR}/qc_${RUN_ACCESSION}_2.fastq.gz"
  local QC="${TMP_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
  local METAPHLAN_LOG="${TMP_DIR}/${RUN_ACCESSION}_metaphlan_log.txt"
  local METAPHLAN_BOWTIE="${TMP_DIR}/${RUN_ACCESSION}_meta.bowtie2out.txt"
  local METAPHLAN_PROFILE="${TMP_DIR}/${RUN_ACCESSION}_metaphlan4.txt"
  local METAPHLAN_COUNTS="${OUTPUT_DIR}/${RUN_ACCESSION}_metaphlan4_counts.txt"
  local MOTUS_PROFILE="${OUTPUT_DIR}/${RUN_ACCESSION}_motus.txt"
  local DADA2_SEQTAB="${OUTPUT_DIR}/seqtab_${RUN_ACCESSION}.rds"
  mkdir -p "$OUTPUT_DIR"

  # Initialize log file
  touch "$LOG_FILE"

  # Validate input files
  if [[ -f "$INPUT_FASTQ" && -f "$PAIRED_FASTQ" ]]; then
    if ! validate_fastq "$INPUT_FASTQ" || ! validate_fastq "$PAIRED_FASTQ"; then
      log_debug "Invalid FASTQ files for $RUN_ACCESSION: $INPUT_FASTQ or $PAIRED_FASTQ"
      append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    fi
    if grep -q "$INPUT_FASTQ" "$LOG_FILE" && grep -q "$PAIRED_FASTQ" "$LOG_FILE"; then
      log_debug "Skipping $RUN_ACCESSION: files already logged"
    else
      append_with_lock "$INPUT_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
      append_with_lock "$PAIRED_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
    fi
  elif [[ -f "$INPUT_FASTQ" ]]; then
    if ! validate_fastq "$INPUT_FASTQ"; then
      log_debug "Invalid FASTQ file for $RUN_ACCESSION: $INPUT_FASTQ"
      append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    fi
    if grep -q "$INPUT_FASTQ" "$LOG_FILE"; then
      log_debug "Skipping $RUN_ACCESSION: file already logged"
    else
      append_with_lock "$INPUT_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
    fi
  else
    log_debug "No valid FASTQ files for $RUN_ACCESSION in $FASTQ_DIR"
    append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  # Determine resources (for logging)
  local FILE_SIZE_GB=$(du -BG "$INPUT_FASTQ" | cut -f1 | sed 's/G//')
  local CPUS=$THREADS_PER_WORKER
  local MEM TIME
  if (( FILE_SIZE_GB <= 5 )); then
    MEM="8G"; TIME="12:00:00"
  elif (( FILE_SIZE_GB <= 20 )); then
    MEM="16G"; TIME="24:00:00"
  else
    MEM="32G"; TIME="48:00:00"
  fi
  log_debug "Processing $RUN_ACCESSION (Type: $SAMPLE_TYPE, Threads: $CPUS)"

  # Workflow
  if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    if [[ -f "$PAIRED_FASTQ" ]]; then
      if [[ ! -f "$QC1" || ! -f "$QC2" ]]; then
        run_command "micromamba run -n \"$DADA2_ENV\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER" \
          "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
      fi
      if validate_fastq "$QC1" && validate_fastq "$QC2" && [[ ! -f "$DADA2_SEQTAB" ]]; then
        run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV\" Rscript run_dada2_partial.R \"$QC1\" \"$QC2\")" \
          "[dada2] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB" 2>/dev/null || log_debug "No seqtab output for $RUN_ACCESSION"
      fi
    else
      if [[ ! -f "$QC" ]]; then
        run_command "micromamba run -n \"$DADA2_ENV\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER" \
          "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
      fi
      if validate_fastq "$QC" && [[ ! -f "$DADA2_SEQTAB" ]]; then
        run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV\" Rscript run_dada2_partial.R \"$QC\")" \
          "[dada2] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB" 2>/dev/null || log_debug "No seqtab output for $RUN_ACCESSION"
      fi
    fi
  elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    if [[ ! -f "$QC" ]]; then
      run_command "micromamba run -n \"$MPA_ENV\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER" \
        "[fastp] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
    fi
    if validate_fastq "$QC" && [[ ! -f "$METAPHLAN_PROFILE" ]]; then
      run_command_with_output "micromamba run -n \"$MPA_ENV\" metaphlan \"$QC\" --input_type fastq --unclassified_estimation --nproc $THREADS_PER_WORKER --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\"" \
        "[metaphlan] Failed for $RUN_ACCESSION" "$METAPHLAN_LOG" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
    fi
    if validate_fastq "$QC" && [[ -f "$METAPHLAN_PROFILE" && ! -f "$METAPHLAN_COUNTS" ]]; then
      convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS" || \
        { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
    fi
    if validate_fastq "$QC" && [[ ! -f "$MOTUS_PROFILE" ]]; then
      run_command "micromamba run -n \"$MOTUS_ENV\" motus profile -s \"$QC\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
        "[motus] Failed for $RUN_ACCESSION" || { append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"; rm -rf "$TMP_DIR"; return 1; }
    fi
  else
    log_debug "Invalid sample type: $SAMPLE_TYPE for $RUN_ACCESSION"
    append_with_lock "$RUN_ACCESSION" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  # Mark as complete only if all steps succeeded
  append_with_lock "${RUN_ACCESSION}:COMPLETE" "$COMPLETED_FILE" "$COMPLETED_LOCK"
  log_debug "Finished processing $RUN_ACCESSION"

  # Clean up temporary directory
  rm -rf "$TMP_DIR"
}

# Final validation and merging
final_validation_and_merge() {
  log_info "Starting final validation and merging"

  # Load expected accessions from input.tsv
  declare -A EXPECTED_SAMPLES
  declare -A BIOPROJECTS
  while IFS="$DELIMITER" read -r bioproject accession sample_type; do
    EXPECTED_SAMPLES["$accession"]="$bioproject $sample_type"
    BIOPROJECTS["$bioproject"]=1
  done < <(tail -n +2 "$INPUT_FILE")

  # Find actual processed samples
  declare -A ACTUAL_SAMPLES
  for dir in */; do
    if [[ -d "$dir" && "$dir" != "*/" ]]; then
      for file in "$dir"*_metaphlan4_counts.txt "$dir"*_motus.txt "$dir"/seqtab_*.rds; do
        if [[ -f "$file" ]]; then
          local accession=$(basename "$file" | sed -E 's/^(metaphlan4_counts|motus|seqtab)_([^_.]+).*/\2/')
          ACTUAL_SAMPLES["$accession"]="$dir"
        fi
      done
    fi
  done

  # Identify missing samples
  local missing=()
  for accession in "${!EXPECTED_SAMPLES[@]}"; do
    if [[ -z "${ACTUAL_SAMPLES[$accession]}" || -z "${COMPLETED_SAMPLES[$accession]}" ]]; then
      missing+=("$accession")
      log_debug "Missing or incomplete sample: $accession"
    fi
  done

  # Retry missing samples
  if [[ ${#missing[@]} -gt 0 ]]; then
    log_info "Retrying ${#missing[@]} missing samples"
    printf '%s\n' "${missing[@]}" | parallel --jobs "$NUM_WORKERS" process_sample "${EXPECTED_SAMPLES[{1}]}" {1}
    load_completed  # Reload to include newly completed samples
  else
    log_info "No missing samples to retry"
  fi

  # Merge profiles for each Bioproject
  for bioproject in "${!BIOPROJECTS[@]}"; do
    merge_profiles "$bioproject" "metaphlan" || log_debug "Failed to merge MetaPhlAn profiles for $bioproject"
    merge_profiles "$bioproject" "motus" || log_debug "Failed to merge mOTUs profiles for $bioproject"
    merge_profiles "$bioproject" "dada2" || log_debug "Failed to merge DADA2 sequence tables for $bioproject"
  done

  # Log final status
  log_info "Final validation and merging complete. Expected: ${#EXPECTED_SAMPLES[@]}, Actual: ${#ACTUAL_SAMPLES[@]}, Completed: ${#COMPLETED_SAMPLES[@]}"
}

# Export functions and variables for GNU Parallel
export -f process_sample log_debug log_info append_with_lock validate_fastq run_command run_command_with_output convert_metaphlan_to_counts merge_profiles final_validation_and_merge generate_quality_profiles
export FASTQ_DIR DEBUG_FILE DEBUG_LOCK COMPLETED_FILE COMPLETED_LOCK FAILED_FILE FAILED_LOCK THREADS_PER_WORKER MOTUS_TAX_LEVEL INPUT_FILE COMPLETED_SAMPLES LOG_LEVEL TMP_BASE SLURM_JOB_ID DELIMITER QUALITY_PROFILE_DIR DADA2_ENV MOTUS_ENV MPA_ENV

# Initialize lock files
touch "$DEBUG_LOCK" "$COMPLETED_LOCK" "$FAILED_LOCK"

# Load completed samples once
load_completed

# Generate quality profiles if requested
if [[ "$QUALITY_CHECK" == "true" ]]; then
  log_info "Generating quality profiles for up to 5 samples"
  generate_quality_profiles
  exit 0
fi

# Process input file with GNU Parallel
log_info "Starting initial processing with $NUM_WORKERS workers"
tail -n +2 "$INPUT_FILE" | parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" process_sample {1} {2} {3}

# Run final validation and merging
final_validation_and_merge

log_info "All processing, validation, and merging complete."
echo "All processing, validation, and merging complete."
