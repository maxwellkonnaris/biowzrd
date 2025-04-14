#!/bin/bash
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

# Purpose: SLURM script for processing 16S and metagenomic FASTQ files using fastp, DADA2, MetaPhlAn, and mOTUs.
# Inputs: CSV/TSV with Bioproject, RunAccession, SequencingType, optional Fastq1/Fastq2.
# Dependencies: micromamba, fastp, R (DADA2), MetaPhlAn, mOTUs, GNU parallel.
# Tip: Use --debug with a small dataset to test the pipeline.

########################################
# Default Variables & Directories
########################################

DEFAULT_DIR="fastq_data/fastq_biologicaldata"
LOCK_DIR="$PWD/locks"
mkdir -p "$LOCK_DIR"

touch "$LOCK_DIR/debug.lock" "$LOCK_DIR/completed.lock" "$LOCK_DIR/failed.lock" "$LOCK_DIR/input.lock"
DEBUG_FILE="debug.log"
DEBUG_LOCK="$LOCK_DIR/debug.lock"
COMPLETED_FILE="completed_steps.log"
COMPLETED_LOCK="$LOCK_DIR/completed.lock"
FAILED_FILE="failed.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
INPUT_LOCK="$LOCK_DIR/input.lock"

DEFAULT_WORKERS=4
DEFAULT_MOTUS_TAX_LEVEL="mOTU"
LOG_LEVEL="INFO"
TMP_BASE="$PWD"
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"
QUALITY_PROFILE_DIR="quality_profiles"

DEFAULT_DADA2_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/dada2"
DEFAULT_MOTUS_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/motus"
DEFAULT_MPA_ENV="/storage/work/mak6930/applicationstorage/micromamba/envs/mpa"

INPUT_HEADER=""
declare -a TMP_DIRS

########################################
# Cleanup
########################################
cleanup() {
  for dir in "${TMP_DIRS[@]}"; do
    [[ -d "$dir" ]] && rm -rf "$dir"
  done
  rm -rf "$LOCK_DIR"
}
trap cleanup EXIT

########################################
# Logging Functions
########################################

log_debug() {
  [[ "$LOG_LEVEL" == "DEBUG" ]] || return
  local msg="$1"
  echo "$(date) [DEBUG] $msg" >&2
  {
    flock -x "$DEBUG_LOCK" || true
    echo "$(date) [DEBUG] $msg" >> "$DEBUG_FILE"
  } 2>/dev/null
}

log_info() {
  local msg="$1"
  echo "$(date) [INFO] $msg" | tee -a "$DEBUG_FILE"
}

########################################
# CSV / File Locking Helpers
########################################

append_with_lock() {
  local line="$1"
  local file="$2"
  local lock="$3"
  {
    flock -x "$lock" || true
    echo "$line" >> "$file"
  } 2>/dev/null
}

update_input_csv() {
  local content="$1"
  local file="$INPUT_FILE"
  local lock="$INPUT_LOCK"
  {
    flock -x "$lock" || true
    printf "%s" "$content" > "$file"
  } 2>/dev/null
}

########################################
# Command-Line Argument Parsing
########################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_FILE="$2"; shift 2 ;;
    -d) FASTQ_DIR="$2"; shift 2 ;;
    -w) NUM_WORKERS="$2"; shift 2 ;;
    -k) MOTUS_TAX_LEVEL="$2"; shift 2 ;;
    -q) QUALITY_CHECK="true"; shift ;;
    -a) DADA2_ENV="$2"; shift 2 ;;
    -m) MOTUS_ENV="$2"; shift 2 ;;
    -p) MPA_ENV="$2"; shift 2 ;;
    --debug) LOG_LEVEL="DEBUG"; shift ;;
    -*) echo "Invalid option: $1" >&2; exit 1 ;;
    *) echo "Unexpected argument: $1" >&2; exit 1 ;;
  esac
done

FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_DIR}"
NUM_WORKERS="${NUM_WORKERS:-$DEFAULT_WORKERS}"
MOTUS_TAX_LEVEL="${MOTUS_TAX_LEVEL:-$DEFAULT_MOTUS_TAX_LEVEL}"
DADA2_ENV="${DADA2_ENV:-$DEFAULT_DADA2_ENV}"
MOTUS_ENV="${MOTUS_ENV:-$DEFAULT_MOTUS_ENV}"
MPA_ENV="${MPA_ENV:-$DEFAULT_MPA_ENV}"

if command -v realpath &>/dev/null; then
  FASTQ_DIR="$(realpath "$FASTQ_DIR")"
fi

########################################
# Check micromamba, environments
########################################
if ! command -v micromamba &>/dev/null; then
  echo "Error: micromamba not found"
  exit 1
fi

check_env() {
  local env="$1"
  local env_name="$2"
  local env_basename

  if [[ -d "$env" ]]; then
    env_basename=$(basename "$env")
    eval "${env_name}_NAME=$env_basename"
  elif micromamba env list | grep -qE "^[[:space:]]*$env[[:space:]]"; then
    eval "${env_name}_NAME=$env"
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

########################################
# Validate input file presence
########################################
if [[ -z "$INPUT_FILE" ]]; then
  echo "Usage: $0 -i <input.tsv|input.csv|input.txt> [-d <fastq_directory>] [-w <num_workers>] [-k <motus_tax_level>] [-q] ..."
  echo "Tip: Use --debug with a small dataset to test the pipeline."
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: Input file $INPUT_FILE does not exist."
  exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ directory $FASTQ_DIR does not exist."
  exit 1
fi

if [[ ! -d "$TMP_BASE" || ! -w "$TMP_BASE" ]]; then
  echo "Error: Temporary directory base $TMP_BASE does not exist or is not writable."
  exit 1
fi

########################################
# Detect input file format
########################################
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

  INPUT_HEADER="$first_line"

  if [[ "$first_line" != "$EXPECTED_HEADER" && "$first_line" != "${EXPECTED_HEADER},Fastq1,Fastq2"* ]]; then
    echo "Error: Input file $file has incorrect header. Expected: $EXPECTED_HEADER or extended with Fastq1,Fastq2"
    exit 1
  fi
}

validate_input_file "$INPUT_FILE"

########################################
# Check RDP DB if 16S present
########################################
check_rdp_database() {
  local has_16s=false
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    if [[ "$sample_type" == "16S" ]]; then
      has_16s=true
      break
    fi
  done < <(tail -n +2 "$INPUT_FILE")
  
  if [[ "$has_16s" == "true" ]]; then
    if [[ ! -f "$RDP_DATABASE" || ! -r "$RDP_DATABASE" ]]; then
      echo "Error: RDP database $RDP_DATABASE not found or not readable."
      exit 1
    fi
    log_info "RDP database $RDP_DATABASE found for 16S samples."
  fi
}
check_rdp_database

########################################
# SLURM / CPU / Worker Setup
########################################
SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
if [[ -z "$NUM_WORKERS" || "$NUM_WORKERS" -le 0 ]]; then
  NUM_WORKERS=$DEFAULT_WORKERS
fi
THREADS_PER_WORKER=$(( SLURM_CPUS_PER_TASK / NUM_WORKERS ))
(( THREADS_PER_WORKER < 1 )) && THREADS_PER_WORKER=1

log_info "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK, NUM_WORKERS=$NUM_WORKERS, THREADS_PER_WORKER=$THREADS_PER_WORKER"
log_info "Running with $NUM_WORKERS workers, $THREADS_PER_WORKER threads each"

########################################
# FASTQ Validation
########################################
validate_fastq() {
  local file="$1"
  local min_size=100
  if [[ ! -f "$file" ]]; then
    echo "File does not exist: $file"
    return 1
  fi

  local size
  size=$(stat -c %s "$file" 2>/dev/null || wc -c < "$file")
  if (( size < min_size )); then
    echo "File too small: $file ($size bytes)"
    return 1
  fi

  if [[ "$file" == *.gz ]]; then
    gzip -t "$file" 2>/dev/null
    if [[ $? -ne 0 ]]; then
      echo "Invalid gzip: $file"
      return 1
    fi
  fi
  return 0
}

########################################
# Command Runners (with retries)
########################################
run_command() {
  local cmd="$1"
  local msg="$2"
  local input_file="$3"
  local output_file="$4"

  local max_attempts=3
  local attempt=1
  local wait_time=2
  local timeout=7200
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command (attempt $attempt): $cmd"
    error_output=$(timeout $timeout bash -c "$cmd" 2>&1)
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
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
  local start_time=$(date +%s)

  while (( attempt <= max_attempts )); do
    log_debug "Running command with output (attempt $attempt): $cmd > $output_file"
    error_output=$(timeout $timeout bash -c "$cmd > \"$output_file\" 2>&1")
    if [[ $? -eq 0 ]]; then
      local end_time=$(date +%s)
      log_info "Command succeeded: [${msg/Failed/Succeeded}] Input: $input_file, Output: $output_file, Time: $((end_time - start_time))s"
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

########################################
# Update Input CSV with FASTQ paths
########################################
update_input_with_fastq_paths() {
  log_info "Checking if input file needs FASTQ path updates"

  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")

  local needs_update=false
  if [[ "$header" == *",Fastq1,Fastq2"* ]]; then
    local all_valid=true
    while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
      if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
        continue
      fi
      if [[ -z "$fastq1" && -z "$fastq2" ]]; then
        all_valid=false
        break
      fi
      if [[ -n "$fastq1" && ! -f "$fastq1" ]]; then
        all_valid=false
        break
      fi
      if [[ -n "$fastq2" && ! -f "$fastq2" ]]; then
        all_valid=false
        break
      fi
    done < "$INPUT_FILE"
    if [[ "$all_valid" == "true" ]]; then
      log_info "FASTQ columns exist and are valid for all samples."
      rm -f "$temp_file"
      return
    else
      needs_update=true
    fi
  else
    header="${header},Fastq1,Fastq2"
    needs_update=true
  fi

  if [[ "$needs_update" == "true" ]]; then
    echo "$header" > "$temp_file"

    declare -A fastq1_map
    declare -A fastq2_map
    while IFS= read -r file; do
      filename=$(basename "$file")
      # Paired-end pattern (_1, _2)
      if [[ "$filename" =~ ^([A-Za-z0-9._-]+)_([12])\.fastq\.gz$ ]]; then
        local accession="${BASH_REMATCH[1]}"
        local pair="${BASH_REMATCH[2]}"
        if [[ "$pair" == "1" ]]; then
          fastq1_map["$accession"]="$file"
        else
          fastq2_map["$accession"]="$file"
        fi
      # Single-end pattern (<accession>.fastq.gz or <accession>_3.fastq.gz)
      elif [[ "$filename" =~ ^([A-Za-z0-9._-]+)(\.fastq\.gz|_3\.fastq\.gz)$ ]]; then
        local accession="${BASH_REMATCH[1]}"
        # Only assign to fastq1 if not already assigned (prioritize paired-end)
        if [[ -z "${fastq1_map[$accession]}" && -z "${fastq2_map[$accession]}" ]]; then
          fastq1_map["$accession"]="$file"
        fi
      fi
    done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq.gz")

    while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
      if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
        continue
      fi
      local new_fastq1="$fastq1"
      local new_fastq2="$fastq2"
      if [[ -z "$fastq1" || ! -f "$fastq1" ]]; then
        new_fastq1="${fastq1_map[$accession]}"
      fi
      if [[ -z "$fastq2" || ! -f "$fastq2" ]]; then
        new_fastq2="${fastq2_map[$accession]}"
      fi
      echo "${bioproject}${DELIMITER}${accession}${DELIMITER}${sample_type}${DELIMITER}${new_fastq1}${DELIMITER}${new_fastq2}${DELIMITER}${rest}" >> "$temp_file"
    done < "$INPUT_FILE"

    local missing_count=0
    while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
      if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
        continue
      fi
      if [[ -z "$fastq1" && -z "$fastq2" ]]; then
        ((missing_count++))
      fi
    done < "$temp_file"

    if [[ $missing_count -gt 0 ]]; then
      log_info "Found $missing_count samples still missing FASTQ files"
      rm -f "$temp_file"
      exit 1
    fi

    update_input_csv "$(cat "$temp_file")"
    rm -f "$temp_file"
    log_info "Updated input file with FASTQ paths"
  fi
}

########################################
# Initialize/Check Checkpoint Columns
########################################
initialize_checkpoints() {
  log_info "Initializing checkpoint columns"
  local temp_file
  temp_file=$(mktemp)
  local header
  header=$(head -n 1 "$INPUT_FILE")

  if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
    if [[ "$header" != *",Fastp,Dada2,Motus,Metaphlan,Completed"* ]]; then
      header="${header},Fastp,Dada2,Motus,Metaphlan,Completed"
    fi
  else
    if [[ "$header" != *",Completed"* ]]; then
      header="${header},Completed"
    fi
  fi
  echo "$header" > "$temp_file"

  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
      continue
    fi
    if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
      if [[ "$header" == *",Fastp,Dada2,Motus,Metaphlan,Completed"* ]]; then
        echo "${bioproject}${DELIMITER}${accession}${DELIMITER}${sample_type}${DELIMITER}${fastq1}${DELIMITER}${fastq2}${DELIMITER}${rest}" >> "$temp_file"
      else
        echo "${bioproject}${DELIMITER}${accession}${DELIMITER}${sample_type}${DELIMITER}${fastq1}${DELIMITER}${fastq2}${DELIMITER}0${DELIMITER}0${DELIMITER}0${DELIMITER}0${DELIMITER}0" >> "$temp_file"
      fi
    else
      local completed="0"
      completed=$(echo "$rest" | awk -F"$DELIMITER" '{print $NF}')
      echo "${bioproject}${DELIMITER}${accession}${DELIMITER}${sample_type}${DELIMITER}${fastq1}${DELIMITER}${fastq2}${DELIMITER}${completed}" >> "$temp_file"
    fi
  done < "$INPUT_FILE"

  update_input_csv "$(cat "$temp_file")"
  rm -f "$temp_file"
  log_info "Checkpoint columns initialized"
}

update_checkpoint() {
  local accession="$1"
  local step="$2"
  local value="$3"
  local temp_file
  temp_file=$(mktemp)

  local header
  header=$(head -n 1 "$INPUT_FILE")
  echo "$header" > "$temp_file"

  local field_index
  case "$step" in
    Fastp) field_index=4 ;;
    Dada2) field_index=5 ;;
    Motus) field_index=6 ;;
    Metaphlan) field_index=7 ;;
    Completed) field_index=8 ;;
  esac
  if [[ "$LOG_LEVEL" != "DEBUG" && "$step" == "Completed" ]]; then
    field_index=6
  fi

  local updated=false
  while IFS="$DELIMITER" read -r bioproject acc sample_type fastq1 fastq2 rest; do
    if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
      continue
    fi
    if [[ "$acc" == "$accession" ]]; then
      local fields=("$bioproject" "$acc" "$sample_type" "$fastq1" "$fastq2" $rest)
      fields[$field_index]="$value"
      echo "$(IFS="$DELIMITER"; echo "${fields[*]}")" >> "$temp_file"
      updated=true
    else
      echo "${bioproject}${DELIMITER}${acc}${DELIMITER}${sample_type}${DELIMITER}${fastq1}${DELIMITER}${fastq2}${DELIMITER}${rest}" >> "$temp_file"
    fi
  done < "$INPUT_FILE"

  if [[ "$updated" == "true" ]]; then
    update_input_csv "$(cat "$temp_file")"
    log_debug "Updated checkpoint for $accession: $step=$value"
  fi
  rm -f "$temp_file"
}

########################################
# Generate Quality Profiles if Requested
########################################
generate_quality_profiles() {
  local max_samples=5
  local count=0
  mkdir -p "$QUALITY_PROFILE_DIR"
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2; do
    if [[ "$sample_type" == "16S" && $count -lt $max_samples ]]; then
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
        fi
      fi
    fi
  done < <(tail -n +2 "$INPUT_FILE")
}

########################################
# Convert Metaphlan profile to counts
########################################
convert_metaphlan_to_counts() {
  local metaphlan_log="$1"
  local metaphlan_profile="$2"
  local output_file="$3"

  if [[ ! -f "$metaphlan_profile" ]]; then
    echo "Missing $metaphlan_profile"
    return 1
  fi

  local mapped_reads
  mapped_reads=$(grep "Total number of reads mapped" "$metaphlan_log" 2>/dev/null | awk '{print $6}' | sed 's/(.*//')
  if [[ -z "$mapped_reads" ]]; then
    mapped_reads=100000
  fi

  awk -v mapped="$mapped_reads" '
    BEGIN { FS="\t"; OFS="\t"; print "#clade_name\trelative_abundance\tread_count" }
    /^#/ { next }
    { count = ($2 * mapped / 100); print $1, $2, count }
  ' "$metaphlan_profile" > "$output_file"

  log_debug "Converted $metaphlan_profile to counts in $output_file"
}

########################################
# Merge Profiles
########################################
merge_profiles() {
  local bioproject="$1"
  local tool="$2"
  local output_dir="$bioproject"
  local merged_file="${output_dir}/${bioproject}_${tool}_merged.txt"

  declare -A expected_accessions
  while IFS="$DELIMITER" read -r proj accession sample_type fastq1 fastq2 rest; do
    if [[ "$proj" == "$bioproject" ]]; then
      expected_accessions["$accession"]=1
    fi
  done < <(tail -n +2 "$INPUT_FILE")

  local all_complete=true
  for accession in "${!expected_accessions[@]}"; do
    local completed
    completed=$(awk -F"$DELIMITER" -v acc="$accession" '$2 == acc {print $NF}' "$INPUT_FILE")
    if [[ "$completed" != "1" ]]; then
      all_complete=false
      break
    fi
  done

  if [[ "$all_complete" == "true" ]]; then
    local profile_files=()
    case "$tool" in
      metaphlan)   profile_files=("$output_dir"/*_metaphlan4_counts.txt) ;;
      motus)       profile_files=("$output_dir"/*_motus.txt) ;;
      dada2)       profile_files=("$output_dir"/seqtab_*.rds) ;;
    esac

    local valid_files=()
    for f in "${profile_files[@]}"; do
      [[ -f "$f" ]] && valid_files+=("$f")
    done

    if [[ "$tool" == "metaphlan" && ${#valid_files[@]} -gt 1 ]]; then
      run_command "micromamba run -n \"$MPA_ENV_NAME\" merge_metaphlan_tables.py ${valid_files[*]} > \"$merged_file\"" \
        "[metaphlan merge] for $bioproject" "${valid_files[*]}" "$merged_file"
    elif [[ "$tool" == "metaphlan" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"
    elif [[ "$tool" == "motus" && ${#valid_files[@]} -gt 1 ]]; then
      local input_list
      input_list=$(printf "%s," "${valid_files[@]}" | sed 's/,$//')
      run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus merge -i \"$input_list\" -o \"$merged_file\"" \
        "[motus merge] for $bioproject" "$input_list" "$merged_file"
    elif [[ "$tool" == "motus" && ${#valid_files[@]} -eq 1 ]]; then
      cp "${valid_files[0]}" "$merged_file"
    elif [[ "$tool" == "dada2" && ${#valid_files[@]} -gt 0 ]]; then
      if [[ ! -f "merge_dada2.R" ]]; then
        log_info "Error: merge_dada2.R not found for DADA2 merging"
        return 1
      fi
      run_command "micromamba run -n \"$DADA2_ENV_NAME\" Rscript merge_dada2.R $bioproject ${valid_files[*]}" \
        "[dada2 merge] for $bioproject" "${valid_files[*]}" ""
    fi
  fi
}

########################################
# Global Array for Completed Samples
########################################
declare -A COMPLETED_SAMPLES
load_completed() {
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    local completed_val
    completed_val=$(echo "$rest" | awk -F"$DELIMITER" '{print $NF}')
    if [[ "$completed_val" == "1" ]]; then
      COMPLETED_SAMPLES["$accession"]=1
    fi
  done < "$INPUT_FILE"
  log_debug "Loaded ${#COMPLETED_SAMPLES[@]} completed samples"
}

########################################
# Process a Single Sample
########################################
process_sample() {
  local BIOPROJECT="$1"
  local RUN_ACCESSION="$2"
  local SAMPLE_TYPE="$3"
  local INPUT_FASTQ="$4"
  local PAIRED_FASTQ="$5"

  log_debug "Processing $RUN_ACCESSION ($SAMPLE_TYPE)"
  
  if [[ "$LOG_LEVEL" == "DEBUG" && -f "$INPUT_FASTQ" ]]; then
    local read_count
    read_count=$(zcat "$INPUT_FASTQ" 2>/dev/null | wc -l)
    read_count=$((read_count/4))
    log_info "$RUN_ACCESSION: Input read count = $read_count"
  fi

  if [[ -n "${COMPLETED_SAMPLES[$RUN_ACCESSION]}" ]]; then
    log_debug "Skipping $RUN_ACCESSION: already completed"
    return
  fi

  local OUTPUT_DIR="${BIOPROJECT}"
  mkdir -p "$OUTPUT_DIR"
  local TMP_DIR
  TMP_DIR=$(mktemp -d -p "$TMP_BASE" "process_${SLURM_JOB_ID}_${RUN_ACCESSION}_XXXXXX")
  TMP_DIRS+=("$TMP_DIR")

  local QC1="${TMP_DIR}/qc_${RUN_ACCESSION}_1.fastq.gz"
  local QC2="${TMP_DIR}/qc_${RUN_ACCESSION}_2.fastq.gz"
  local QC="${TMP_DIR}/qc_${RUN_ACCESSION}.fastq.gz"
  local METAPHLAN_LOG="${TMP_DIR}/${RUN_ACCESSION}_metaphlan_log.txt"
  local METAPHLAN_BOWTIE="${TMP_DIR}/${RUN_ACCESSION}_meta.bowtie2out.txt"
  local METAPHLAN_PROFILE="${TMP_DIR}/${RUN_ACCESSION}_metaphlan4.txt"
  local METAPHLAN_COUNTS="${OUTPUT_DIR}/${RUN_ACCESSION}_metaphlan4_counts.txt"
  local MOTUS_PROFILE="${OUTPUT_DIR}/${RUN_ACCESSION}_motus.txt"
  local DADA2_SEQTAB="${OUTPUT_DIR}/seqtab_${RUN_ACCESSION}.rds"

  if [[ -f "$INPUT_FASTQ" ]]; then
    local reason1
    reason1=$(validate_fastq "$INPUT_FASTQ")
    if [[ $? -ne 0 ]]; then
      append_with_lock "$RUN_ACCESSION:$reason1" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    fi
  else
    local reason="Missing FASTQ: $INPUT_FASTQ"
    append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  if [[ -n "$PAIRED_FASTQ" && -f "$PAIRED_FASTQ" ]]; then
    local reason2
    reason2=$(validate_fastq "$PAIRED_FASTQ")
    if [[ $? -ne 0 ]]; then
      append_with_lock "$RUN_ACCESSION:$reason2" "$FAILED_FILE" "$FAILED_LOCK"
      rm -rf "$TMP_DIR"
      return 1
    fi
  fi

  # =============== 16S Workflow ===============
  if [[ "$SAMPLE_TYPE" == "16S" ]]; then
    if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
      local fastp_status
      fastp_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $4}' "$INPUT_FILE")
      if [[ "$fastp_status" == "1" ]]; then
        log_debug "Skipping fastp for $RUN_ACCESSION: checkpoint set"
      else
        if [[ -n "$PAIRED_FASTQ" ]]; then
          run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
            "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2"
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
        else
          run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
            "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC"
          update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
        fi
      fi
    else
      if [[ -n "$PAIRED_FASTQ" ]]; then
        run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2"
      else
        run_command "micromamba run -n \"$DADA2_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC"
      fi
    fi

    if [[ -n "$PAIRED_FASTQ" ]]; then
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$QC1\" \"$QC2\")" \
        "[dada2] Process for $RUN_ACCESSION" "$QC1,$QC2" "$DADA2_SEQTAB"
      if [[ -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
      else
        log_info "Error: DADA2 output $TMP_DIR/seqtab_${RUN_ACCESSION}.rds not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 output missing" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    else
      run_command "(cd \"$TMP_DIR\" && micromamba run -n \"$DADA2_ENV_NAME\" Rscript run_dada2_partial.R \"$QC\")" \
        "[dada2] Process for $RUN_ACCESSION" "$QC" "$DADA2_SEQTAB"
      if [[ -f "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" ]]; then
        mv "$TMP_DIR/seqtab_${RUN_ACCESSION}.rds" "$DADA2_SEQTAB"
      else
        log_info "Error: DADA2 output $TMP_DIR/seqtab_${RUN_ACCESSION}.rds not found for $RUN_ACCESSION"
        append_with_lock "$RUN_ACCESSION:DADA2 output missing" "$FAILED_FILE" "$FAILED_LOCK"
        rm -rf "$TMP_DIR"
        return 1
      fi
    fi

  # =============== Metagenomic Workflow ===============
  elif [[ "$SAMPLE_TYPE" == "meta" ]]; then
    if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
      local fastp_status
      fastp_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $4}' "$INPUT_FILE")
      if [[ "$fastp_status" != "1" ]]; then
        if [[ -n "$PAIRED_FASTQ" ]]; then
          run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
            "[fastp] Process,
            "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2"
        else
          run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
            "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC"
        fi
        update_checkpoint "$RUN_ACCESSION" "Fastp" "1"
      fi
    else
      if [[ -n "$PAIRED_FASTQ" ]]; then
        run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -I \"$PAIRED_FASTQ\" -o \"$QC1\" -O \"$QC2\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ,$PAIRED_FASTQ" "$QC1,$QC2"
      else
        run_command "micromamba run -n \"$MPA_ENV_NAME\" fastp -i \"$INPUT_FASTQ\" -o \"$QC\" -w $THREADS_PER_WORKER --disable_quality_filtering --empty_output_tabular" \
          "[fastp] Process for $RUN_ACCESSION" "$INPUT_FASTQ" "$QC"
      fi
    fi

    local FASTQ_TO_USE
    if [[ -n "$PAIRED_FASTQ" ]]; then
      FASTQ_TO_USE="$QC1,$QC2"
    else
      FASTQ_TO_USE="$QC"
    fi

    if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
      local metaphlan_status
      metaphlan_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $7}' "$INPUT_FILE")
      if [[ "$metaphlan_status" != "1" ]]; then
        run_command "micromamba run -n \"$MPA_ENV_NAME\" metaphlan \"$FASTQ_TO_USE\" --input_type fastq --unclassified_estimation --nproc $THREADS_PER_WORKER --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\" 2> \"$METAPHLAN_LOG\"" \
          "[metaphlan] Process for $RUN_ACCESSION" "$FASTQ_TO_USE" "$METAPHLAN_PROFILE"
        if [[ -f "$METAPHLAN_PROFILE" ]]; then
          convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS"
        fi
        update_checkpoint "$RUN_ACCESSION" "Metaphlan" "1"
      fi
    else
      if [[ ! -f "$METAPHLAN_PROFILE" ]]; then
        run_command "micromamba run -n \"$MPA_ENV_NAME\" metaphlan \"$FASTQ_TO_USE\" --input_type fastq --unclassified_estimation --nproc $THREADS_PER_WORKER --bowtie2out \"$METAPHLAN_BOWTIE\" -o \"$METAPHLAN_PROFILE\" 2> \"$METAPHLAN_LOG\"" \
          "[metaphlan] Process for $RUN_ACCESSION" "$FASTQ_TO_USE" "$METAPHLAN_PROFILE"
        if [[ -f "$METAPHLAN_PROFILE" ]]; then
          convert_metaphlan_to_counts "$METAPHLAN_LOG" "$METAPHLAN_PROFILE" "$METAPHLAN_COUNTS"
        fi
      fi
    fi

    if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
      local motus_status
      motus_status=$(awk -F"$DELIMITER" -v acc="$RUN_ACCESSION" '$2 == acc {print $6}' "$INPUT_FILE")
      if [[ "$motus_status" != "1" ]]; then
        if [[ -n "$PAIRED_FASTQ" ]]; then
          run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -f \"$QC1\" -r \"$QC2\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
            "[motus] Process for $RUN_ACCESSION" "$QC1,$QC2" "$MOTUS_PROFILE"
        else
          run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -s \"$QC\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
            "[motus] Process for $RUN_ACCESSION" "$QC" "$MOTUS_PROFILE"
        fi
        update_checkpoint "$RUN_ACCESSION" "Motus" "1"
      fi
    else
      if [[ ! -f "$MOTUS_PROFILE" ]]; then
        if [[ -n "$PAIRED_FASTQ" ]]; then
          run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -f \"$QC1\" -r \"$QC2\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
            "[motus] Process for $RUN_ACCESSION" "$QC1,$QC2" "$MOTUS_PROFILE"
        else
          run_command "micromamba run -n \"$MOTUS_ENV_NAME\" motus profile -s \"$QC\" -o \"$MOTUS_PROFILE\" -t $THREADS_PER_WORKER -c -k $MOTUS_TAX_LEVEL" \
            "[motus] Process for $RUN_ACCESSION" "$QC" "$MOTUS_PROFILE"
        fi
      fi
    fi
  else
    local reason="Invalid sample type: $SAMPLE_TYPE"
    append_with_lock "$RUN_ACCESSION:$reason" "$FAILED_FILE" "$FAILED_LOCK"
    rm -rf "$TMP_DIR"
    return 1
  fi

  update_checkpoint "$RUN_ACCESSION" "Completed" "1"
  append_with_lock "${RUN_ACCESSION}:COMPLETE" "$COMPLETED_FILE" "$COMPLETED_LOCK"
  log_info "Finished $RUN_ACCESSION"
  rm -rf "$TMP_DIR"
}

########################################
# Final Merge & Validation
########################################
final_validation_and_merge() {
  log_info "Starting final validation and merging"

  declare -A EXPECTED_SAMPLES
  declare -A BIOPROJECTS

  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
      continue
    fi
    EXPECTED_SAMPLES["$accession"]="$bioproject $sample_type"
    BIOPROJECTS["$bioproject"]=1
  done < <(tail -n +2 "$INPUT_FILE")

  for bioproject in "${!BIOPROJECTS[@]}"; do
    merge_profiles "$bioproject" "metaphlan"
    merge_profiles "$bioproject" "motus"
    merge_profiles "$bioproject" "dada2"
  done

  log_info "Final validation and merging complete."
  log_info "All processing complete."
}

########################################
# Main Workflow
########################################

update_input_with_fastq_paths
initialize_checkpoints
load_completed

if [[ "$QUALITY_CHECK" == "true" ]]; then
  log_info "Generating quality profiles for up to 5 samples"
  generate_quality_profiles
  exit 0
fi

cat << 'EOF' > /tmp/process_sample_wrapper.sh
#!/bin/bash
$(declare -f process_sample)
$(declare -f log_debug)
$(declare -f log_info)
$(declare -f append_with_lock)
$(declare -f update_checkpoint)
$(declare -f validate_fastq)
$(declare -f run_command)
$(declare -f run_command_with_output)
$(declare -f convert_metaphlan_to_counts)
BIOPROJECT="$1"
RUN_ACCESSION="$2"
SAMPLE_TYPE="$3"
INPUT_FASTQ="$4"
PAIRED_FASTQ="$5"
DEBUG_FILE="$6"
DEBUG_LOCK="$7"
COMPLETED_FILE="$8"
COMPLETED_LOCK="$9"
FAILED_FILE="${10}"
FAILED_LOCK="${11}"
THREADS_PER_WORKER="${12}"
MOTUS_TAX_LEVEL="${13}"
INPUT_FILE="${14}"
LOG_LEVEL="${15}"
TMP_BASE="${16}"
SLURM_JOB_ID="${17}"
DELIMITER="${18}"
DADA2_ENV_NAME="${19}"
MOTUS_ENV_NAME="${20}"
MPA_ENV_NAME="${21}"
process_sample "$BIOPROJECT" "$RUN_ACCESSION" "$SAMPLE_TYPE" "$INPUT_FASTQ" "$PAIRED_FASTQ"
EOF
chmod +x /tmp/process_sample_wrapper.sh

log_info "Starting parallel processing with $NUM_WORKERS workers"

if [[ "$LOG_LEVEL" == "DEBUG" ]]; then
  log_info "DEBUG mode: Processing only one 16S and one meta sample sequentially"
  local processed_16s=0
  local processed_meta=0
  while IFS="$DELIMITER" read -r bioproject accession sample_type fastq1 fastq2 rest; do
    if [[ -z "$bioproject" || "$bioproject" == "Bioproject" ]]; then
      continue
    fi
    if [[ "$sample_type" == "16S" && $processed_16s -eq 0 ]]; then
      /tmp/process_sample_wrapper.sh "$bioproject" "$accession" "$sample_type" "$fastq1" "$fastq2" \
        "$DEBUG_FILE" "$DEBUG_LOCK" "$COMPLETED_FILE" "$COMPLETED_LOCK" "$FAILED_FILE" "$FAILED_LOCK" \
        "$THREADS_PER_WORKER" "$MOTUS_TAX_LEVEL" "$INPUT_FILE" "$LOG_LEVEL" "$TMP_BASE" "$SLURM_JOB_ID" \
        "$DELIMITER" "$DADA2_ENV_NAME" "$MOTUS_ENV_NAME" "$MPA_ENV_NAME"
      ((processed_16s++))
    elif [[ "$sample_type" == "meta" && $processed_meta -eq 0 ]]; then
      /tmp/process_sample_wrapper.sh "$bioproject" "$accession" "$sample_type" "$fastq1" "$fastq2" \
        "$DEBUG_FILE" "$DEBUG_LOCK" "$COMPLETED_FILE" "$COMPLETED_LOCK" "$FAILED_FILE" "$FAILED_LOCK" \
        "$THREADS_PER_WORKER" "$MOTUS_TAX_LEVEL" "$INPUT_FILE" "$LOG_LEVEL" "$TMP_BASE" "$SLURM_JOB_ID" \
        "$DELIMITER" "$DADA2_ENV_NAME" "$MOTUS_ENV_NAME" "$MPA_ENV_NAME"
      ((processed_meta++))
    fi
    if [[ $processed_16s -eq 1 && $processed_meta -eq 1 ]]; then
      break
    fi
  done < "$INPUT_FILE"
else
  tail -n +2 "$INPUT_FILE" | \
    parallel --colsep "$DELIMITER" --jobs "$NUM_WORKERS" --halt now,fail=1 \
      /tmp/process_sample_wrapper.sh {1} {2} {3} {4} {5} \
      "$DEBUG_FILE" "$DEBUG_LOCK" "$COMPLETED_FILE" "$COMPLETED_LOCK" "$FAILED_FILE" "$FAILED_LOCK" \
      "$THREADS_PER_WORKER" "$MOTUS_TAX_LEVEL" "$INPUT_FILE" "$LOG_LEVEL" "$TMP_BASE" "$SLURM_JOB_ID" \
      "$DELIMITER" "$DADA2_ENV_NAME" "$MOTUS_ENV_NAME" "$MPA_ENV_NAME"
fi

rm -f /tmp/process_sample_wrapper.sh

final_validation_and_merge

log_info "All processing, validation, and merging complete."
echo "All processing, validation, and merging complete."
