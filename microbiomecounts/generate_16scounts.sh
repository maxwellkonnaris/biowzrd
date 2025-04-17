#!/bin/bash
#SBATCH --job-name=counts_16s
#SBATCH --output=slurm_16s-%j.out
#SBATCH --error=slurm_16s-%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --account=one
#SBATCH --mail-user=mak6930@psu.edu
#SBATCH --open-mode=append
#--------------------------------------------------
# Purpose : SLURM pipeline for 16S‐rRNA FASTQ files
# Tooling : micromamba + R (DADA2) + GNU parallel
#--------------------------------------------------

set -euxo pipefail

########################################
# Default variables & directories
########################################
DEFAULT_DIR="fastq_data/fastq_biologicaldata"
OUTPUT_BASE="$PWD/dada2_16s"
LOCK_DIR="$PWD/locks"
FAILED_FILE="failed_16s.log"
FAILED_LOCK="$LOCK_DIR/failed.lock"
INPUT_LOCK="$LOCK_DIR/input.lock"
DEFAULT_WORKERS=4
LOG_LEVEL="INFO"
RDP_DATABASE="rdp_19_toGenus_trainset.fa.gz"
REPAIR_FASTQS=1

mkdir -p "$LOCK_DIR" "$OUTPUT_BASE" || {
  echo "Error: cannot create required directories"; exit 1; }

touch "$FAILED_LOCK" "$INPUT_LOCK"

########################################
# Cleanup
########################################
cleanup() { rm -rf "$LOCK_DIR"; }
trap cleanup EXIT

########################################
# Logging helpers
########################################
log_debug() {
  if [[ $LOG_LEVEL == DEBUG ]]; then
    echo "$(date) [DEBUG] $1" >&2
  fi
}
log_info()  { echo  "$(date) [INFO]  $1"; }

########################################
# File‐locking helpers
########################################
append_with_lock() {          # $1 line  $2 file  $3 lock
  { flock -x 200; echo "$1" >> "$2"; } 200>>"$3"; }

update_input_csv() {          # $1 new‐file‐contents
  { flock -x 200; printf '%s\n' "$1" > "$INPUT_FILE"; } 200>"$INPUT_LOCK"; }

########################################
# Argument parsing
########################################
PROCESS_SAMPLE=0
while [[ $# -gt 0 ]]; do
  case $1 in
    -i) INPUT_FILE=$2; shift 2 ;;
    -d) FASTQ_DIR=$2; shift 2 ;;
    -w) NUM_WORKERS=$2; shift 2 ;;
    --debug) LOG_LEVEL=DEBUG; shift ;;
    --no-repair-fastqs) REPAIR_FASTQS=0; shift ;;
    --process-sample) PROCESS_SAMPLE=1; PROCESS_ARGS=("$@"); break ;;
    *) echo "Invalid option: $1"; exit 1;;
  esac
done

FASTQ_DIR=${FASTQ_DIR:-$DEFAULT_DIR}
NUM_WORKERS=${NUM_WORKERS:-$DEFAULT_WORKERS}

command -v realpath &>/dev/null && FASTQ_DIR=$(realpath "$FASTQ_DIR")

########################################
# Conda / micromamba environment
########################################
command -v micromamba >/dev/null || { echo "micromamba not in PATH"; exit 1; }

# Assume DADA2 environment is activated
DADA2_ENV_NAME=$(micromamba env list | grep '\*' | awk '{print $1}')
if [[ -z "$DADA2_ENV_NAME" ]]; then
  echo "Error: No active micromamba environment detected"
  exit 1
fi
log_info "Using activated DADA2 environment: $DADA2_ENV_NAME"

########################################
# Generic cmd‐runner under the env
########################################
run_command() {               # $1 cmd  $2 description
  local cmd="$1" desc="$2"
  log_info "$desc"
  micromamba run -n "$DADA2_ENV_NAME" -- bash -c "$cmd"
}

########################################
# Validate CSV header, set delimiter
########################################
validate_input_file() {
  local first
  first=$(head -n1 "$1")
  if [[ $first == "Bioproject,RunAccession,SequencingType"* ]]; then
    DELIM=','; EXPECT="Bioproject,RunAccession,SequencingType"
  elif [[ $first == $'Bioproject\tRunAccession\tSequencingType'* ]]; then
    DELIM=$'\t'; EXPECT=$'Bioproject\tRunAccession\tSequencingType'
  else
    echo "Error: malformed header"; exit 1;
  fi
  [[ $first == "$EXPECT"* ]] || { echo "Error: header mismatch"; exit 1; }
}
validate_input_file "$INPUT_FILE"
export DELIM

########################################
# Check reference database
########################################
[[ -r $RDP_DATABASE ]] || { echo "Missing RDP database: $RDP_DATABASE"; exit 1; }
log_info "RDP DB ok: $RDP_DATABASE"

########################################
# FASTQ helpers
########################################
generate_fastq_checksums() {
  log_info "Generating initial checksums for FASTQ files"
  find "$FASTQ_DIR" -type f -name "*.fastq.gz" | sort | \
    parallel --jobs "$NUM_WORKERS" 'b3sum "{}"' \
    > "$FASTQ_DIR/checksums.b3"
  log_info "Checksums saved to $FASTQ_DIR/checksums.b3"
}

validate_fastq() {
  local fq=$1
  [[ -f $fq ]] || return 1
  gzip -t "$fq" 2>/dev/null || return 1

  zcat "$fq" | head -n 3 | awk '
    NR==1 && $0 !~ /^@/ { print "[FAIL] Line 1 does not start with @"; exit 1 }
    NR==3 && $0 !~ /^\+/ { print "[FAIL] Line 3 does not start with +"; exit 1 }
    END { print "[PASS] FASTQ header structure OK" }
  ' || return 1

  return 0
}



repair_fastq_if_needed() {
  local src_dir="$FASTQ_DIR"
  local repaired_dir="${src_dir%/}_repaired"
  local checksum_file="$src_dir/checksums.b3"
  local failed_log="$src_dir/failed_checksums.txt"
  local old_pwd="$PWD"

  mkdir -p "$repaired_dir"
  : > "$failed_log"

  log_info "Verifying FASTQ checksums using b3sum in parallel …"

  cd "$src_dir"

  # Generate checksums if not present
  [[ -f "$checksum_file" ]] || generate_fastq_checksums

  # 1) Build “to-check” list (hash <TAB> file)
  build_validated_set
  printf '%s\n' "${VALIDATED_ACCESSIONS[@]}" | sed 's/$/.fastq.gz/' > skip_files.txt
  grep -vFf skip_files.txt "$checksum_file" > checksums_to_check.b3
  rm -f skip_files.txt

  # 2) Parallel b3sum verification
  awk '{print $1 "\t" $2}' checksums_to_check.b3 | \
    parallel --colsep '\t' \
             --jobs "$NUM_WORKERS" --eta \
             'echo -e "{1}  {2}" | b3sum -c - --quiet || echo "{2}"' \
    > "$failed_log"

  rm -f checksums_to_check.b3

  local n_failed
  n_failed=$(wc -l < "$failed_log" | tr -d ' ')
  log_info "$n_failed files failed checksum"

  # 3) Repair corrupted FASTQs
  if (( n_failed > 0 )); then
    log_info "Repairing corrupted FASTQ files …"
    parallel --jobs "$NUM_WORKERS" --eta '
      fq="{}"
      fq_base="${fq##*/}"
      fq_out="'"$repaired_dir"'/$fq_base"
      echo "[REPAIR] $fq_base" >&2
      zcat "$fq" | sed "s/\r$//" | pigz -p 4 > "$fq_out"
    ' < "$failed_log"
  fi

  # 4) Link untouched files
  log_info "Linking valid FASTQ files …"
  find "$src_dir" -type f -name "*.fastq.gz" -print0 | \
    while IFS= read -r -d '' fq; do
      fq_base="${fq##*/}"
      [[ -e "$repaired_dir/$fq_base" ]] || ln -s "$(realpath "$fq")" "$repaired_dir/$fq_base"
    done

  # 5) Regenerate BLAKE3 checksums
  log_info "Writing new checksums …"
  (cd "$repaired_dir" && b3sum *.fastq.gz > checksums.b3)

  # 6) Update FASTQ_DIR
  FASTQ_DIR="$repaired_dir"
  export FASTQ_DIR
  cd "$old_pwd" >/dev/null
  log_info "Repair complete → FASTQ_DIR updated to $FASTQ_DIR"
}

is_biological_fastq() {
  local fq=$1
  local lengths
  lengths=$(zcat "$fq" | awk 'NR % 4 == 2 {print length}' | head -n 20)
  local uniq
  uniq=$(printf '%s\n' "$lengths" | sort -u | wc -l)
  local maxlen
  maxlen=$(printf '%s\n' "$lengths" | awk 'max < $1 {max = $1} END {print max}')
  # biological if:   ≥2 different lengths   OR   any read > 30 bp
  (( uniq > 1 )) || (( maxlen > 30 ))
}

mark_failure() {              # $1 code  $2 accession
  append_with_lock "$2:$1" "$FAILED_FILE" "$FAILED_LOCK"
}

filter_completed_samples() {
  log_info "Filtering completed samples from input"

  local header nfields comp_index
  header=$(head -n1 "$INPUT_FILE")
  IFS="$DELIM" read -ra fields <<< "$header"

  # Find the index of the Completed column
  for i in "${!fields[@]}"; do
    if [[ "${fields[$i]}" == "Completed" ]]; then
      comp_index=$((i+1))
      break
    fi
  done

  if [[ -z "$comp_index" ]]; then
    log_info "No 'Completed' column found — keeping all rows"
    return
  fi

  {
    echo "$header"
    tail -n +2 "$INPUT_FILE" | awk -F"$DELIM" -v idx="$comp_index" '$idx != 1'
  } > tmp_input.$$

  mv tmp_input.$$ "$INPUT_FILE"
  log_info "Filtered completed rows. Remaining: $(($(wc -l < "$INPUT_FILE") - 1)) samples"
}

build_validated_set() {
  VALIDATED_ACCESSIONS=()
  local val_index
  local header
  header=$(head -n1 "$INPUT_FILE")
  IFS="$DELIM" read -ra fields <<< "$header"
  for i in "${!fields[@]}"; do
    if [[ "${fields[$i]}" == "Validated" ]]; then
      val_index=$((i+1))
      break
    fi
  done
  [[ -z "$val_index" ]] && return

  mapfile -t VALIDATED_ACCESSIONS < <(
    tail -n +2 "$INPUT_FILE" | awk -F"$DELIM" -v idx="$val_index" '$idx == 1 {print $2}'
  )
}

########################################
# Build FASTQ index
########################################
find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | sort > fastq_index_16s.txt

########################################
# Check for unmatched accessions (optional diagnostic)
########################################
log_info "Checking for unmatched accessions between FASTQ files and input CSV"

# 1. Extract accessions from CSV (column 2, skip header)
cut -d',' -f2 "$INPUT_FILE" | grep -v RunAccession | sort > 16saccessions_from_csv.txt

# 2. Extract accessions from FASTQ index by filename
awk -F/ '{fn = $NF; gsub(/(_[12])?\.fastq(\.gz)?$/, "", fn); print fn}' fastq_index_16s.txt | sort > 16saccessions_from_index.txt

# 3. Find accessions in CSV not found in FASTQ filenames
comm -23 16saccessions_from_csv.txt 16saccessions_from_index.txt > unmatched_16saccessions.txt

# 4. Report count
missing_count=$(wc -l < unmatched_16saccessions.txt)
log_info "$missing_count accessions in the input CSV were not found in any FASTQ file"
if (( missing_count > 0 )); then
  log_info "Example unmatched accessions:"
  head -n 10 unmatched_16saccessions.txt >&2
fi

rm -rf 16saccessions_from_csv.txt 16saccessions_from_index.txt

########################################
# Update input with FASTQ paths
########################################
update_input_with_fastq_paths() {
  log_info "Updating FASTQ paths"

  # ---------------------------------------------------------------------
  # 1) Fix header (strip CR, add missing columns if needed)
  # ---------------------------------------------------------------------
  local header new_header
  header=$(head -n1 "$INPUT_FILE" | tr -d '\r')

  if [[ "$header" == *",Fastq1,Fastq2,Validated"* ]]; then
    new_header="$header"
  else
    new_header="$header,Fastq1,Fastq2,Validated"
  fi
  printf '%s\n' "$new_header" > tmp.$$

  # ---------------------------------------------------------------------
  # 2) Build map: accession → "fq1 fq2 ..."
  # ---------------------------------------------------------------------
  declare -A fastq_map
  while IFS= read -r path; do
    fn=$(basename "$path")
    acc=$(echo "$fn" | sed -E 's/(_[1-4])?\.fastq(\.gz)?$//' | xargs)
    fastq_map["$acc"]+="${path} "
  done < fastq_index_16s.txt

  # ---------------------------------------------------------------------
  # 3) Loop over CSV rows (avoid subshell to preserve array)
  # ---------------------------------------------------------------------
  tail -n +2 "$INPUT_FILE" > lines_to_process.$$
  while IFS="$DELIM" read -r biop acc st f1 f2 val rest || [[ -n "$biop" ]]; do
    acc=$(echo "$acc" | tr -d '\r' | xargs)
    st=$(echo "$st" | tr -d '\r' | xargs)

    if [[ "$st" != "16S" ]]; then
      printf '%s\n' "$biop$DELIM$acc$DELIM$st$DELIM$f1$DELIM$f2$DELIM$val$DELIM$rest" >> tmp.$$
      continue
    fi

    if [[ -n "${acc:-}" ]]; then
      echo "[DEBUG] Acc: $acc → Files: ${fastq_map[$acc]:-}" >&2
      if [[ -z "${fastq_map[$acc]:-}" ]]; then
        echo "[WARN] No FASTQ files found for accession $acc" >&2
      fi
    fi

    IFS=' ' read -r -a files <<< "${fastq_map[$acc]:-}"

    valid=()
    for fq in "${files[@]}"; do
      [[ -z "$fq" ]] && continue
      echo "[DEBUG] Evaluating $fq for $acc" >&2
    
      if ! [[ -f "$fq" ]]; then
        echo "[FAIL] $fq does not exist" >&2
        mark_failure MISSING "$acc"
        continue
      fi
    
      if ! gzip -t "$fq" 2>/dev/null; then
        echo "[FAIL] $fq failed gzip test" >&2
        mark_failure CORRUPT "$acc"
        continue
      fi
    
      if ! zcat "$fq" | head -n 3 | awk '
          NR==1 && $0 !~ /^@/ { print "[FAIL] Line 1 not @"; exit 1 }
          NR==3 && $0 !~ /^\+/ { print "[FAIL] Line 3 not +"; exit 1 }'; then
        echo "[FAIL] $fq failed FASTQ structure" >&2
        mark_failure CORRUPT "$acc"
        continue
      fi
    
      lengths=$(zcat "$fq" | awk 'NR % 4 == 2 {print length}' | head -n 20)
      uniq=$(printf '%s\n' "$lengths" | sort -u | wc -l)
      maxlen=$(printf '%s\n' "$lengths" | awk 'max < $1 {max = $1} END {print max}')
    
      echo "[DEBUG] $fq → uniq=$uniq, maxlen=$maxlen" >&2
    
      if (( uniq > 1 || maxlen > 30 )); then
        echo "[PASS] $fq is biological" >&2
        valid+=("$fq")
      else
        echo "[FAIL] $fq is technical" >&2
        mark_failure TECHNICAL "$acc"
      fi
    done
    
    echo "[RESULT] $acc → ${#valid[@]} valid files: ${valid[*]}" >&2


    case ${#valid[@]} in
      0) new1="MISSING"; new2="MISSING"; val=0; mark_failure NO_VALID "$acc" ;;
      1) new1="${valid[0]}"; new2="";    val=1 ;;
      *) new1="${valid[0]}"; new2="${valid[1]}"; val=1 ;;
    esac

    printf '%s\n' "$biop$DELIM$acc$DELIM$st$DELIM$new1$DELIM$new2$DELIM$val$DELIM$rest" >> tmp.$$
  done < lines_to_process.$$
  rm lines_to_process.$$

  # ---------------------------------------------------------------------
  # 4) Overwrite input under flock
  # ---------------------------------------------------------------------
  mv tmp.$$ "$INPUT_FILE"
  rm tmp.$$
  log_info "FASTQ path update done"
}

########################################
# Checkpoint columns
########################################
initialize_checkpoints() {
  log_info "Ensuring checkpoint columns"
  local hdr
  hdr=$(head -n1 "$INPUT_FILE")
  [[ $hdr == *",Validated,Dada2,Completed"* ]] || hdr="$hdr,Validated,Dada2,Completed"
  {
    echo "$hdr"
    tail -n +2 "$INPUT_FILE" | while IFS="$DELIM" read -r b a st f1 f2 v d2 c rest; do
      [[ -z $b || $b == "Bioproject" ]] && continue
      echo "$b$DELIM$a$DELIM$st$DELIM$f1$DELIM$f2$DELIM${v:-0}$DELIM${d2:-0}$DELIM${c:-0}"
    done
  } > tmp.$$
  update_input_csv "$(cat tmp.$$)"
  rm tmp.$$
  log_info "Checkpoint columns ok"
}

update_checkpoint() {         # $1 accession  $2 field  $3 value
  local accession=$1 field=$2 value=$3 idx
  case $field in
    Validated) idx=5 ;; Dada2) idx=6 ;; Completed) idx=7 ;; *) return ;;
  esac
  {
    read -r hdr
    echo "$hdr"
    while IFS="$DELIM" read -r b a st f1 f2 v d2 c rest; do
      [[ -z $b || $b == "Bioproject" ]] && continue
      if [[ $a == "$accession" ]]; then
        fields=("$b" "$a" "$st" "$f1" "$f2" "$v" "$d2" "$c")
        fields[$idx]=$value
        printf '%s' "${fields[0]}"
        for i in {1..7}; do printf '%s%s' "$DELIM" "${fields[$i]}"; done
        printf '\n'
      else
        echo "$b$DELIM$a$DELIM$st$DELIM$f1$DELIM$f2$DELIM$v$DELIM$d2$DELIM$c"
      fi
    done
  } < "$INPUT_FILE" > tmp.$$
  update_input_csv "$(cat tmp.$$)"
  rm tmp.$$
}

########################################
# process_sample (worker)
########################################
process_sample() {
  local BIOP=$1 ACC=$2 TYPE=$3 FQ1=$4 FQ2=$5 VAL=$6 D2=$7 COMP=$8
  [[ $TYPE != 16S ]] && { log_debug "Skip $ACC (not 16S)"; return; }
  [[ "$VAL" == 0 ]] && { log_debug "Skip $ACC (invalid FASTQ)"; return; }
  [[ -z $FQ1 && -z $FQ2 ]] && { mark_failure "NO_FASTQ" "$ACC"; return; }
  [[ "$COMP" == 1 ]] && { log_debug "Already done $ACC"; return; }

  local outdir="$OUTPUT_BASE/$BIOP"
  mkdir -p "$outdir"

  local seqtab="${outdir}/asv_${ACC}.rds"
  if [[ "$D2" != 1 ]]; then
    if [[ -n $FQ2 ]]; then
      cmd="Rscript \"$PWD/run_dada2_partial.R\" \"$FQ1\" \"$FQ2\" \"$seqtab\""
    else
      cmd="Rscript \"$PWD/run_dada2_partial.R\" \"$FQ1\" \"$seqtab\""
    fi
    run_command "$cmd" "DADA2 on $ACC"
    [[ -f $seqtab ]] || { mark_failure "DADA2_FAIL" "$ACC"; return; }
    update_checkpoint "$ACC" Dada2 1
  fi

  update_checkpoint "$ACC" Completed 1
  log_info "Finished $ACC"
}
export -f process_sample log_debug log_info append_with_lock update_checkpoint run_command validate_fastq is_biological_fastq
export DELIM INPUT_FILE FAILED_FILE FAILED_LOCK INPUT_LOCK OUTPUT_BASE DADA2_ENV_NAME LOG_LEVEL

########################################
# Merging per‐bioproject profiles
########################################
merge_profiles() {
  local biop=$1 odir="$OUTPUT_BASE/$biop"
  local done=1 acc
  while IFS="$DELIM" read -r b a st _ _ _ _ comp; do
    [[ $b == "$biop" && $st == "16S" ]] || continue
    [[ $comp == 1 ]] || { done=0; break; }
  done < <(tail -n +2 "$INPUT_FILE")
  (( done )) || return

  mapfile -t tabs < <(ls "$odir"/asv_*.rds 2>/dev/null)
  [[ ${#tabs[@]} -gt 0 ]] || return
  run_command "Rscript \"$PWD/merge_dada2.R\" $biop ${tabs[*]}" "Merging $biop"
}

final_validation_and_merge() {
  log_info "Running final validation"
  awk -F"$DELIM" 'NR>1 && $3=="16S"{print $1}' "$INPUT_FILE" | sort -u |
    while read -r bp; do merge_profiles "$bp"; done
  log_info "Final validation done"
}

########################################
# Parallel driver
########################################
process_samples() {
  log_info "Launching GNU parallel ($NUM_WORKERS workers)"
  tail -n +2 "$INPUT_FILE" | awk -F"$DELIM" '$3=="16S"' | \
  parallel --colsep "$DELIM" \
           --jobs "$NUM_WORKERS" \
           --halt 2,fail=1 \
           --joblog "$OUTPUT_BASE/parallel_joblog.txt" \
           'process_sample {1} {2} {3} {4} {5} {6} {7} {8}'
  log_info "All samples finished"
}

########################################
# Main
########################################
main() {
  [[ -f "$FASTQ_DIR/checksums.b3" ]] || generate_fastq_checksums
  (( REPAIR_FASTQS )) && repair_fastq_if_needed
  update_input_with_fastq_paths
  initialize_checkpoints
  process_samples
  final_validation_and_merge
  log_info "Pipeline complete."
}

if (( PROCESS_SAMPLE )); then
  process_sample "${PROCESS_ARGS[@]}"
else
  main
fi
