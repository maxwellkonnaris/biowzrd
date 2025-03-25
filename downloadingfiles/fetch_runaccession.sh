#!/bin/bash
#SBATCH --time=00:30:00   # Adjust as needed
#SBATCH --mem=1G          # Adjust memory as needed
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err

set -euo pipefail

# ---[ Basic sanity checks for commands we rely on ]---
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

# ---[ Get the study ID from the command line ]---
STUDY_ID=${1:-}
if [[ -z "$STUDY_ID" ]]; then
    echo "Usage: sbatch $0 <STUDY_ID>"
    exit 1
fi

echo "🔍 Processing study: $STUDY_ID"

# ---[ Define the output file and a separate lock file ]---
OUTPUT_FILE="run_accessions.txt"
LOCK_FILE="run_accessions.lock"

# Make sure OUTPUT_FILE exists so grep won't complain
touch "$OUTPUT_FILE"

# Acquire an exclusive lock (via file descriptor 200)
exec 200>"$LOCK_FILE"
flock -x 200   # Blocks until lock is available

# ---[ Start of Critical Section ]---
append_unique_runs() {
    local runs="$1"
    if [[ -z "$runs" ]]; then
        echo "⚠️ No runs found for $STUDY_ID"
        return 1
    fi

    # Check for duplicates before appending
    while IFS= read -r run; do
        if ! grep -qx "${run}" "$OUTPUT_FILE"; then
            echo "$run" >> "$OUTPUT_FILE"
        fi
    done <<< "$runs"
}

# ---[ Process based on study type ]---
case "$STUDY_ID" in
    PRJNA*|SRP*)
        # NCBI SRA query (using EDirect)
        runs=$(esearch -db sra -query "$STUDY_ID" \
            | efetch -format runinfo \
            | awk -F, 'NR>1 {print $1}' \
            | sort -u)
        append_unique_runs "$runs"
        ;;
    PRJEB*|ERP*|EGAS*)
        # ENA query (via API)
        runs=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY_ID&result=read_run&fields=run_accession&format=tsv" \
            | awk 'NR>1 {print $1}' \
            | sort -u)
        append_unique_runs "$runs"
        ;;
    GSE*)
        # NCBI GEO → SRA conversion (limited, may need manual inspection)
        runs=$(esearch -db gds -query "$STUDY_ID" \
            | elink -target sra \
            | efetch -format runinfo \
            | awk -F, 'NR>1 {print $1}' \
            | sort -u)
        append_unique_runs "$runs"
        ;;
    *)
        echo "❌ Unsupported study type: $STUDY_ID (Skipping)"
        exit 1
        ;;
esac
# ---[ End of Critical Section ]---

echo "✅ Finished processing $STUDY_ID"
