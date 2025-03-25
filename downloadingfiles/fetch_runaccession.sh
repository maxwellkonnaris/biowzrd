#!/bin/bash
#SBATCH --time=06:00:00        # Adjust as needed
#SBATCH --mem=2G               # Adjust as needed
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --job-name=fetch_runs
#SBATCH --output=logs/fetch_runs_all.out
#SBATCH --error=logs/fetch_runs_all.err

set -euo pipefail

# -----------------------------------------------------
# 1) Basic checks for commands we rely on
# -----------------------------------------------------
if ! command -v esearch >/dev/null 2>&1; then
    echo "ERROR: esearch not found in PATH."
    exit 1
fi

if ! command -v efetch >/dev/null 2>&1; then
    echo "ERROR: efetch not found in PATH."
    exit 1
fi

if ! command -v elink >/dev/null 2>&1; then
    echo "ERROR: elink not found in PATH."
    exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl not found in PATH."
    exit 1
fi

# -----------------------------------------------------
# 2) Define input file and outputs
# -----------------------------------------------------
STUDY_FILE="studies.txt"
OUTPUT_FILE="run_accessions.txt"
LOCK_FILE="run_accessions.lock"

# Check that the study file exists
if [[ ! -f "$STUDY_FILE" ]]; then
  echo "ERROR: $STUDY_FILE not found. Please create it with one study ID per line."
  exit 1
fi

# Optional: Uncomment if you want to start fresh every time
# rm -f "$OUTPUT_FILE"

# Ensure the output file exists (so grep won't complain)
touch "$OUTPUT_FILE"

# Acquire an exclusive lock (via file descriptor 200)
exec 200>"$LOCK_FILE"
flock -x 200   # This blocks until the lock is available

# -----------------------------------------------------
# 3) Helper function to append runs uniquely
# -----------------------------------------------------
append_unique_runs() {
    local runs="$1"

    if [[ -z "$runs" ]]; then
        echo "⚠️  No runs found"
        return 1
    fi

    # For each run, append if not already in the file
    while IFS= read -r run; do
        if ! grep -qx "${run}" "$OUTPUT_FILE"; then
            echo "$run" >> "$OUTPUT_FILE"
        fi
    done <<< "$runs"
}

# -----------------------------------------------------
# 4) Main loop: read each study, fetch run accessions
# -----------------------------------------------------
while read -r STUDY_ID; do
    # Skip blank lines
    [[ -z "$STUDY_ID" ]] && continue

    echo "🔍 Processing study: $STUDY_ID"

    # Case: what type of study ID is this?
    case "$STUDY_ID" in
        PRJNA*|SRP*)
            # SRA query via EDirect (CSV)
            runs=$(
                esearch -db sra -query "$STUDY_ID" \
                | efetch -format runinfo \
                | sed 's/\r$//' \
                | awk -F, 'NR>1 {print $1}' \
                | sort -u
            )
            append_unique_runs "$runs"
            ;;
        PRJEB*|ERP*|EGAS*)
            # ENA query (TSV)
            runs=$(
                curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY_ID&result=read_run&fields=run_accession&format=tsv" \
                | sed 's/\r$//' \
                | awk 'NR>1 {print $1}' \
                | sort -u
            )
            append_unique_runs "$runs"
            ;;
        GSE*)
            # NCBI GEO → SRA
            runs=$(
                esearch -db gds -query "$STUDY_ID" \
                | elink -target sra \
                | efetch -format runinfo \
                | sed 's/\r$//' \
                | awk -F, 'NR>1 {print $1}' \
                | sort -u
            )
            append_unique_runs "$runs"
            ;;
        *)
            echo "❌ Unsupported study type: $STUDY_ID (skipping)"
            ;;
    esac

    echo "✅ Done: $STUDY_ID"
    echo
done < "$STUDY_FILE"

# -----------------------------------------------------
# 5) Done!
# -----------------------------------------------------
echo "🏁 All studies processed!"
