#!/bin/bash
# Ensure bash is used
export SHELL=/bin/bash

#SBATCH --job-name=counts
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --account=open
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --mail-type=END,FAIL

# Default directory for FASTQ files
DEFAULT_DIR="fastq_data/fastq_biologicaldata"
LOCK_DIR="$PWD/locks"
mkdir -p "$LOCK_DIR"
DEBUG_FILE="debug.log"
DEBUG_LOCK="$LOCK_DIR/debug.lock"
COMPLETED_FILE="completed_steps.log"
COMPLETED_LOCK="$LOCK_DIR/completed.lock"
FAILED_FILE="failed.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
DEFAULT_WORKERS=4
DEFAULT_MOTUS_TAX_LEVEL="mOTU"
LOG_LEVEL="DEBUG"
TMP_BASE="$PWD"
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"
QUALITY_PROFILE_DIR="quality_profiles"
DEFAULT_DADA2_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/dada2"
DEFAULT_MOTUS_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/motus"
DEFAULT_MPA_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa"

# Cleanup function
cleanup() {
  rm -rf "$LOCK_DIR"
}
trap cleanup EXIT

# Logging functions
log_debug() {
  [[ "$LOG_LEVEL" == "DEBUG" ]] || return
  local msg="$1"
  echo "$(date) [DEBUG] $msg" >&2
  flock -x "$DEBUG_LOCK" -c "echo \"$(date) [DEBUG] $msg\" >> \"$DEBUG_FILE\"" 2>/dev/null || \
    echo "[WARN] Could not log to $DEBUG_FILE: $msg" >&2
}

log_info() {
  local msg="$1"
  echo "$(date) [INFO] $msg" | tee -a "$DEBUG_FILE"
}

# Append to file with locking
append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  flock -x "$lock" -c "echo \"$line\" >> \"$file\"" 2>/dev/null || \
    echo "[ERROR] Could not append '$line' to $file" >&2
}

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

# Normalize FASTQ_DIR path
FASTQ_DIR=$(echo "$FASTQ_DIR" | sed 's|//|/|g')

# Validate micromamba
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found"
  exit 1
fi

# Check environment
check_env() {
  local env="$1"
  local env_name="$2"
  local env_basename

  if [[ -d "$env" ]]; then
    env_basename=$(basename "$env")
    eval "${env_name}_NAME=$env_basename"
    return 0
  elif micromamba env list | grep -qE "^[[:space:]]*$env[[:space:]]"; then
    eval "${env_name}_NAME=$env"
    return 0
  else
    echo "Error: micromamba environment $env_name ($env) not found"
    exit 1
  fi
}

check_env "$DADA2_ENV" "dada2"
check_env "$MOTUS_ENV" "motus"
check_env "$MPA_ENV" "mpa"

DADA2_ENV_NAME="${dada2_NAME}"
MOTUS_ENV_NAME="${motus_NAME}"
MPA_ENV_NAME="${mpa_NAME}"

# Validate input file
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv|input.csv|input.txt> [-d <fastq_directory>] [-w <num_workers>] [-k <motus_tax_level>] [-q] [-a <dada2_env>] [-m <motus_env>] [-p <mpa_env>]"
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file $INPUT_FILE does not exist."
  exit 1
fi

# Validate temporary directory
if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
  echo "Error: Temporary directory base $TMP_BASE does not exist or is not writable."
  exit 1
fi

# Detect input file format
validate_input_file() {
  local file="$1"
  local first_line
  first_line=$(head -n 1 "$file")
  
  if [[ "$first_line" == *"Bioproject,RunAccession,SequencingType"* ]]; then
    DELIMITER=','
    EXPECTED_HEADER="Bioproject,RunAccession,SequencingType"
  elif [[ "$first_line" == $'Bioproject\tRunAccession\tSequencingType'* ]]; then
    DELIMITER=$'\t'
    EXPECTED_HEADER="Bioproject\tRunAccession\tSequencingType"
  else
    echo "Error: Input file $file has unrecognized format or header."
    exit 1
  fi

  if [[ "$first_line" != "$EXPECTED_HEADER" ]]; then
    echo "Error: Input file $file has incorrect header. Expected: $EXPECTED_HEADER"
    exit 1
  fi
}
validate_input_file "$INPUT_FILE"

# Check RDP database
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
      echo "Error: RDP database $RDP_DATABASE not found."
      exit 1
    fi
    log_info "RDP database $RDP_DATABASE found for 16S samples."
  fi
}
check_rdp_database

# Validate FASTQ directory
if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

# Set SLURM defaults
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
NUM_WORKERS=${NUM_WORKERS:-$DEFAULT_WORKERS}

# Calculate threads per worker
if [[ -z "$NUM_WORKERS" || "$NUM_WORKERS" -le 0 ]]; then
  NUM_WORKERS=4
fi
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
if [[ $THREADS_PER_WORKER -lt 1 ]]; then
  THREADS_PER_WORKER=1
  NUM_WORKERS=$SLURM_CPUS_PER_TASK
fi

# Log configuration
log_info "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK, NUM_WORKERS=$NUM_WORKERS, THREADS_PER_WORKER=$THREADS_PER_WORKER"
log_info "Running with $NUM_WORKERS workers, $THREADS_PER_WORKER threads each (total $SLURM_CPUS_PER_TASK CPUs)"

# Validate FASTQ file
validate_fastq() {
  local file="$1"
  local min_size=100
  if [[ ! -f "$file" ]]; then
    log_debug "[validate_fastq] File $file does not exist"
    echo "File does not exist: $file"
    return 1
  fi
  local size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    log_debug "[validate_fastq] File $file too small: $size bytes"
    echo "File too small: $file ($size bytes)"
    return 1
  fi
  if [[ "$file" == *.gz && ! $(gzip -t "$file" 2>/dev/null && echo "valid") ]]; then
    log_debug "[validate_fastq] File $file is not a valid gzip"
    echo "Invalid gzip: $file"
    return 1
  fi
  return 0
}

run_command() {
  local cmd="$1"
  local msg="$2"
  local input_file="$3"
  local output_file="$4"
  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=7200
  local error_output
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command (attempt $attempt): $cmd"
    error_output=$(timeout $timeout bash -c "$cmd" 2>&1)
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [$(echo "$msg" | sed 's/Failed/Succeeded/')] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

run_command_with_output() {
  local cmd="$1"
  local msg="$2"
  local output_file="$3"
  local input_file="$4"
  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=7200
  local error_output
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command with output (attempt $attempt): $cmd > $output_file"
    error_output=$(timeout $timeout bash -c "$cmd > \"$output_file\" 2>&1")
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [$(echo "$msg" | sed 's/Failed/Succeeded/')] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
      return 0
    else
      log_debug "$msg (attempt $attempt/$max_attempts), Input: $input_file, Output: $output_file, Error: $error_output"
      if (( attempt == max_attempts )); then
        log_info "$msg failed after $max_attempts attempts, Input: $input_file, Output: $output_file, Error: $error_output"
        echo "$msg: $error_output"
        return 1
      fi
      sleep $(( wait_time * attempt ))
      (( attempt++ ))
    fi
  done
}

generate_quality_profiles() {
  local max_samples=5
  local count=0
  mkdir -p "$QUALITY_PROFILE_DIR"
  while IFS="$DELIMITER" read -r bioproject accession sample_type; do
    if [[ "$sample_type" == "16S" && $count -lt $max_samples ]]; then
      local fastq1="${FASTQ_DIR}/${accession}_1.fastq.gz"
      local fastq2="${FASTQ_DIR}/${accession}_2.fastq.gz"
      local output_pdf="${QUALITY_PROFILE_DIR}/${accession}_quality.pdf"
      if [[ -f "$fastq1" && ! -f "$output_pdf" ]]; then
        if [[ -f "$fastq2" ]]; then
          run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript -e 'library(dada2); pdf(\"$output_pdf\"); plotQualityProfile(c(\"$fastq1\", \"$fastq2\")); dev.off()'" \
            "[quality profile] Generate for $accession" "$fastq1,$fastq2" "$output_pdf"
        else
          run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript -e 'library(dada2); pdf(\"$output_pdf\"); plotQualityProfile(\"$fastq1\"); dev.off()'" \
            "[quality profile] Generate for $accession" "$fastq1" "$output_pdf"
        fi
        if [[ $? -eq 0 ]]; then
          log_info "Generated quality profile for $accession in $output_pdf"
          (( count++ ))
        else
          log_debug "Failed to generate quality profile for $accession"
        fi
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
    echo "Missing files: $metaphlan_log or $metaphlan_profile"
    return 1
  fi
  local mapped_reads=$(grep "Total number of reads mapped" "$metaphlan_log" | awk '{print $6}' | sed 's/(.*//')
  if [[ -z "$mapped_reads" ]]; then
    log_debug "[convert_metaphlan_to_counts] Could not find mapped reads in $metaphlan_log"
    echo "Could not find mapped reads in $metaphlan_log"
    return 1
  fi
  awk -v mapped="$mapped_reads" '
    BEGIN { FS="\t"; OFS="\t"; print "#clade_name\trelative_abundance\tread_count" }
    /^#/ { next }
    { count = ($2 * mapped / 100); print $1, $2, count }
  ' "$metaphlan_profile" > "$output_file"
  log_debug "[convert_metaphlan_to_counts] Converted $metaphlan_profile to counts in $output_file"
  return 0
}

# Merge profiles
merge_profiles() {
  local bioproject="$1"
  local tool="$2"
  local output_dir="$bioproject"
  local merged_file="${output_dir}/${bioproject}_${tool}_merged.txt"
  local profile_files

  declare -A expected_accessions
  while IFS="$DELIMITER" read -r proj accession sample_type; do
    if [[ "$proj" == "$bioproject" ]]; then
      expected_accessions["$accession"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  if [[ "$tool" == "metaphlan" ]]; then
    profile_files=($(ls "$output_dir"/*_metaphlan4_counts.txt 2>/dev/null))
  elif [[ "$tool" == "motus" ]]; then
    profile_files=($(ls "$output_dir"/*_motus.txt 2>/dev/null))
  elif [[ "$tool" == "dada2" ]]; then
    profile_files=($(ls "$output_dir"/seqtab_*.rds 2>/dev/null))
  fi

  local all_complete=true
  for accession in "${!expected_accessions[@]}"; do
    if [[ -z "${COMPLETED_SAMPLES[$accession]}" ]]; then
      log_debug "[merge_profiles] $bioproject ($tool): Accession $accession not processed"
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

  if [[ "$all_complete" == "true" && ${#profile_files[@]} -gt 0 ]]; then
    if [[ "$tool" == "metaphlan" && ${#profile_files[@]} -gt 1 ]]; then
      local input_list="${profile_files[*]}"
      run_command "micromamba run -n \"$MPA_ENV_NAME\" merge_metaphlan_tables.py $input_list > \"$merged_file\"" \
        "[metaphlan merge] Process for $bioproject" "${profile_files[*]}" "$merged_file" || return 1
      log_debug "Merged MetaPhlAn profiles for $bioproject into $merged_file"
    elif [[ "$tool" == "metaphlan" && ${#profile_files[@]} -eq 1 ]]; then
      cp "${profile_files[0]}" "$merged_file"
      log_debug "Copied single metaphlan profile for $bioproject to $merged_file"
    elif [[ "$tool" == "motus" && ${#profile_files[@]} -gt 1 ]]; then
      local input_list=$(printf "%s," "${profile_files[@]}" | sed 's/,$//')
      run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus merge -i \"$input_list\" -o \"$merged_file\"" \
        "[motus merge] Process for $bioproject" "$input_list" "$merged_file" || return 1
      log_debug "Merged mOTUs profiles for $bioproject into $merged_file"
    elif [[ "$tool" == "motus" && ${#profile_files[@]} -eq 1 ]]; then
      cp "${profile_files[0]}" "$merged_file"
      log_debug "Copied single motus profile for $bioproject to $merged_file"
    elif [[ "$tool" == "dada2" ]]; then
      local input_list="${bioproject} ${profile_files[*]}"
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript merge_dada2.R $input_list" \
        "[dada2 merge] Process for $bioproject" "${profile_files[*]}" "" || return 1
      log_debug "Merged DADA2 sequence tables for $bioproject"
    fi
  else
    log_debug "[merge_profiles] $bioproject ($tool): Skipping merge due to incomplete profiles"
    return 1
  fi
  return 0
}

# Load completed samples
declare -A COMPLETED_SAMPLES
load_completed() {
  if [[ -f "$COMPLETED_FILE" ]]; then
    while IFS=":" read -r accession status; do
      if [[ "$status" == "COMPLETE" ]]; then
        COMPLETED_SAMPLES["$accession"]=1
      fi
    done < "$COMPLETED_FILE"
  fi
  log_debug "Loaded ${#COMPLETED_SAMPLES[@]} completed samples"
}

process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SAMPLE_TYPE="$3"

  log_debug "Starting process_sample for $RUN_ACCESSION ($SAMPLE_TYPE)"
  log_info "Processing $RUN_ACCESSION: Input files $INPUT_FASTQ, $PAIRED_FASTQ"
  if [[ -f "$INPUT_FASTQ" ]]; then
    local read_count=$(zcat "$INPUT_FASTQ" 2>/dev/null | echo $((`wc -l`/4)))
    log_info "$RUN_ACCESSION: Input read count = $read_count"
  fi

  if [[ -n "${COMPLETED_SAMPLES[$RUN_ACCESSION]}" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already processed"
    return
  fi

  local INPUT_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_1.fastq.gz"
  local PAIRED_FASTQ="${FASTQ_DIR}/${RUN_ACCESSION}_2.fastq.gz"
  local OUTPUT_DIR="${BIOPROJECT}"
  local LOG_FILE="${OUTPUT_DIR}/processed_files.log"
  local TMP_DIR=$(mktemp -d -p "$TMP_BASE" "process_${SLURM_JOB_ID}_${RUN_ACCESSION}_XXXXXX")
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

  touch "$LOG_FILE"

  # Validate input files
  if [[ -f "$INPUT_FASTQ" && -f "$PAIRED_FASTQ" ]]; then
    local reason1=$(validate_fastq "$INPUT_FASTQ")
    local status1=$?
    local reason2=$(validate_fastq "$PAIRED_FASTQ")
    local status2=$?
    if [[ $status1 -ne 0 || $status2 -ne 0 ]]; then
      local reason="Invalid FASTQ files - $INPUT_FASTQ: $reason1; $PAIRED_FASTQ: $reason2"
      log_debug "$reason"
      append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
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
    local reason=$(validate_fastq "$INPUT_FASTQ")
    if [[ $? -ne 0 ]]; then
      log_debug "Invalid FASTQ file for $RUN_ACCESSION: $reason"
      append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    fi
    if grep -q "$INPUT_FASTQ" "$LOG_FILE"; then
      log_debug "Skipping $RUN_ACCESSION: file already logged"
    else
      append_with_lock "$INPUT_FASTQ" "$LOG_FILE" "$COMPLETED_LOCK"
    fi
  else
    local reason="No valid FASTQ files for $RUN_ACCESSION in $FASTQ_DIR"
    log_debug "$reason"
    append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  local FILE_SIZE_GB=$(du -BG "$INPUT_FASTQ" | cut -f1 | sed 's/G//')
  log_debug "Processing $RUN_ACCESSION (Type: $SAMPLE_TYPE, Size: ${FILE_SIZE_GB}GB, Threads: $THREADS_PER_WORKER)"

  if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    if [[ -f "$PAIRED_FASTQ" ]]; then
      if [[ ! -f "$QC1" || ! -f "$QC2" ]]; then
        local reason=$(run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER -q 10 -u 80 -l 10 --detect_adapter_for_pe --disable_quality_filtering --trim_poly_g --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
      local reason1=$(validate_fastq "$QC1")
      local status1=$?
      local reason2=$(validate_fastq "$QC2")
      local status2=$?
      if [[ $status1 -ne 0 || $status2 -ne 0 ]]; then
        local reason="Invalid QC FASTQ files - $QC1: $reason1; $QC2: $reason2"
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
      log_info "$RUN_ACCESSION: Post-fastp read count = $(zcat \"$QC1\" 2>/dev/null | echo $((`wc -l`/4)))"
      if [[ ! -f "$DADA2_SEQTAB" ]]; then
        local reason=$(run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$QC1\" \"$QC2\")" \
          "[dada2] Process for $RUN_ACCESSION" "$QC1,$QC2" "$DADA2_SEQTAB")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB" 2>/dev/null || {
          local reason="No seqtab output for $RUN_ACCESSION"
          log_debug "$reason"
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        }
      fi
    else
      if [[ ! -f "$QC" ]]; then
        local reason=$(run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER -q 10 -u 80 -l 10 --trim_poly_g --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
      local reason=$(validate_fastq "$QC")
      if [[ $? -ne 0 ]]; then
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
      log_info "$RUN_ACCESSION: Post-fastp read count = $(zcat \"$QC\" 2>/dev/null | echo $((`wc -l`/4)))"
      if [[ ! -f "$DADA2_SEQTAB" ]]; then
        local reason=$(run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$QC\")" \
          "[dada2] Process for $RUN_ACCESSION" "$QC" "$DADA2_SEQTAB")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB" 2>/dev/null || {
          local reason="No seqtab output for $RUN_ACCESSION"
          log_debug "$reason"
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        }
      fi
    fi
  elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    if [[ -f "$PAIRED_FASTQ" ]]; then
      if [[ ! -f "$QC1" || ! -f "$QC2" ]]; then
        local reason=$(run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER -q 10 -u 80 -l 10 --detect_adapter_for_pe --disable_quality_filtering --trim_poly_g --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
      local reason1=$(validate_fastq "$QC1")
      local status1=$?
      local reason2=$(validate_fastq "$QC2")
      local status2=$?
      if [[ $status1 -ne 0 || $status2 -ne 0 ]]; then
        local reason="Invalid QC FASTQ files - $QC1: $reason1; $QC2: $reason2"
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
      log_info "$RUN_ACCESSION: Post-fastp read count = $(zcat \"$QC1\" 2>/dev/null | echo $((`wc -l`/4)))"
      local FASTQ_TO_USE="$QC1,$QC2"
    else
      if [[ ! -f "$QC" ]]; then
        local reason=$(run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER -q 10 -u 80 -l 10 --trim_poly_g --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC")
        if [[ $? -ne 0 ]]; then
          append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
          rm -rf "$TMP_DIR"
          return 1
        fi
      fi
      local reason=$(validate_fastq "$QC")
      if [[ $? -ne 0 ]]; then
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
      log_info "$RUN_ACCESSION: Post-fastp read count = $(zcat \"$QC\" 2>/dev/null | echo $((`wc -l`/4)))"
      local FASTQ_TO_USE="$QC"
    fi
    if [[ ! -f "$METAPHLAN_PROFILE" ]]; then
      rm -f "$METAPHLAN_BOWTIE"
      local reason=$(run_command_with_output "micromamba run -n \"$MPA_ENV_NAME\" metaphlan \"$FASTQ_TO_USE\" --input_type fastq --unclassified_estimation --nproc $THREADS_PER_WORKER --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\"" \
        "[metaphlan] Process for $RUN_ACCESSION" "$METAPHLAN_LOG" "$FASTQ_TO_USE")
      if [[ $? -ne 0 ]]; then
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    fi
    if [[ -f "$METAPHLAN_PROFILE" && ! -f "$METAPHLAN_COUNTS" ]]; then
      local reason=$(convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS")
      if [[ $? -ne 0 ]]; then
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    fi
    if [[ ! -f "$MOTUS_PROFILE" ]]; then
      if [[ -f "$PAIRED_FASTQ" ]]; then
        local reason=$(run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -f \"$QC1\" -r \"$QC2\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$QC1,$QC2" "$MOTUS_PROFILE")
      else
        local reason=$(run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -s \"$QC\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
          "[motus] Process for $RUN_ACCESSION" "$QC" "$MOTUS_PROFILE")
      fi
      if [[ $? -ne 0 ]]; then
        append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    fi
  else
    local reason="Invalid sample type: $SAMPLE_TYPE"
    log_debug "$reason for $RUN_ACCESSION"
    append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  append_with_lock "${RUN_ACCESSION}:COMPLETE" "$COMPLETED_FILE" "$COMPLETED_LOCK"
  log_info "Finished processing $RUN_ACCESSION"
  rm -rf "$TMP_DIR"
}

# Final validation and merging
final_validation_and_merge() {
  log_info "Starting final validation and merging"

  declare -A EXPECTED_SAMPLES
  declare -A BIOPROJECTS
  while IFS="$DELIMITER" read -r bioproject accession sample_type; do
    EXPECTED_SAMPLES["$accession"]="$bioproject $sample_type"
    BIOPROJECTS["$bioproject"]=1
  done < <(tail -n +2 "$INPUT_FILE")

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

  local missing=()
  for accession in "${!EXPECTED_SAMPLES[@]}"; do
    if [[ -z "${ACTUAL_SAMPLES[$accession]}" || -z "${COMPLETED_SAMPLES[$accession]}" ]]; then
      missing+=("$accession")
      log_debug "Missing or incomplete sample: $accession"
    fi
  done

  if [[ ${#missing[@]} -gt 0 ]]; then
    log_info "Retrying ${#missing[@]} missing samples"
    printf '%s\n' "${missing[@]}" | parallel --jobs "$NUM_WORKERS" --halt now,fail=1 process_sample "${EXPECTED_SAMPLES[{1}]}" {1}
    load_completed
  else
    log_info "No missing samples to retry"
  fi

  for bioproject in "${!BIOPROJECTS[@]}"; do
    merge_profiles "$bioproject" "metaphlan" || log_debug "Failed to merge MetaPhlAn profiles for $bioproject"
    merge_profiles "$bioproject" "motus" || log_debug "Failed to merge mOTUs profiles for $bioproject"
    merge_profiles "$bioproject" "dada2" || log_debug "Failed to merge DADA2 sequence tables for $bioproject"
  done

  log_info "Final validation and merging complete. Expected: ${#EXPECTED_SAMPLES[@]}, Actual: ${#ACTUAL_SAMPLES[@]}, Completed: ${#COMPLETED_SAMPLES[@]}"
}

# Export functions and variables
export -f process_sample log_debug log_info append_with_lock validate_fastq run_command run_command_with_output convert_metaphlan_to_counts merge_profiles final_validation_and_merge generate_quality_profiles
export FASTQ_DIR DEBUG_FILE DEBUG_LOCK COMPLETED_FILE COMPLETED_LOCK FAILED_FILE FAILED_LOCK THREADS_PER_WORKER MOTUS_TAX_LEVEL INPUT_FILE COMPLETED_SAMPLES LOG_LEVEL TMP_BASE SLURM_JOB_ID DELIMITER QUALITY_PROFILE_DIR DADA2_ENV_NAME MOTUS_ENV_NAME MPA_ENV_NAME SLURM_CPUS_PER_TASK NUM_WORKERS

# Initialize lock files
touch "$DEBUG_LOCK" "$COMPLETED_LOCK" "$FAILED_LOCK"

# Load completed samples
load_completed

# Generate quality profiles if requested
if [[ "$QUALITY_CHECK" == "true" ]]; then
  log_info "Generating quality profiles for up to 5 samples"
  generate_quality_profiles
  exit 0
fi

# Create wrapper script for GNU Parallel
cat << 'EOF' > /tmp/process_sample_wrapper.sh
#!/bin/bash
BIOPROJECT="$1"
RUN_ACCESSION="$2"
SAMPLE_TYPE="$3"
export FASTQ_DIR="$4"
export DEBUG_FILE="$5"
export DEBUG_LOCK="$6"
export COMPLETED_FILE="$7"
export COMPLETED_LOCK="$8"
export FAILED_FILE="$9"
export FAILED_LOCK="${10}"
export THREADS_PER_WORKER="${11}"
export MOTUS_TAX_LEVEL="${12}"
export INPUT_FILE="${13}"
export LOG_LEVEL="${14}"
export TMP_BASE="${15}"
export SLURM_JOB_ID="${16}"
export DELIMITER="${17}"
export DADA2_ENV_NAME="${18}"
export MOTUS_ENV_NAME="${19}"
export MPA_ENV_NAME="${20}"
export -f process_sample log_debug log_info append_with_lock validate_fastq run_command run_command_with_output convert_metaphlan_to_counts
process_sample "$BIOPROJECT" "$RUN_ACCESSION" "$SAMPLE_TYPE"
EOF
chmod +x /tmp/process_sample_wrapper.sh

# Process input file
log_info "Starting initial processing with $NUM_WORKERS workers"
log_debug "Parallel command: tail -n +2 \"$INPUT_FILE\" | parallel --colsep \"$DELIMITER\" --jobs $NUM_WORKERS --halt now,fail=1 /tmp/process_sample_wrapper.sh {1} {2} {3} \"$FASTQ_DIR\" \"$DEBUG_FILE\" \"$DEBUG_LOCK\" \"$COMPLETED_FILE\" \"$COMPLETED_LOCK\" \"$FAILED_FILE\" \"$FAILED_LOCK\" \"$THREADS_PER_WORKER\" \"$MOTUS_TAX_LEVEL\" \"$INPUT_FILE\" \"$LOG_LEVEL\" \"$TMP_BASE\" \"$SLURM_JOB_ID\" \"$DELIMITER\" \"$DADA2_ENV_NAME\" \"$MOTUS_ENV_NAME\" \"$MPA_ENV_NAME\""
tail -n +2 "$INPUT_FILE" | parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" --halt now,fail=1 /tmp/process_sample_wrapper.sh {1} {2} {3} "$FASTQ_DIR" "$DEBUG_FILE" "$DEBUG_LOCK" "$COMPLETED_FILE" "$COMPLETED_LOCK" "$FAILED_FILE" "$FAILED_LOCK" "$THREADS_PER_WORKER" "$MOTUS_TAX_LEVEL" "$INPUT_FILE" "$LOG_LEVEL" "$TMP_BASE" "$SLURM_JOB_ID" "$DELIMITER" "$DADA2_ENV_NAME" "$MOTUS_ENV_NAME" "$MPA_ENV_NAME"

# Clean up wrapper
rm -f /tmp/process_sample_wrapper.sh

# Final validation and merging
final_validation_and_merge

log_info "All processing, validation, and merging complete."
echo "All processing, validation, and merging complete."
