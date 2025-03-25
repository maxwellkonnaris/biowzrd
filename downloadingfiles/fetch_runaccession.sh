#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --job-name=fetch_runs
#SBATCH --output=logs/fetch_runs_all.out
#SBATCH --error=logs/fetch_runs_all.err

# -------------------------------------------
# (1) Basic checks for commands
# -------------------------------------------
for cmd in esearch efetch elink curl; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found in PATH."
        exit 1
    fi
done

# -------------------------------------------
# (2) Input/output setup
# -------------------------------------------
STUDY_FILE="studies.txt"
OUTPUT_FILE="run_accessions.txt"

if [[ ! -f "$STUDY_FILE" ]]; then
    echo "ERROR: $STUDY_FILE not found. Please create it with one study ID per line."
    exit 1
fi

touch "$OUTPUT_FILE"

# -------------------------------------------
# (3) Process each study
# -------------------------------------------
COUNT=0

while read -r STUDY_ID || [[ -n "$STUDY_ID" ]]; do
    [[ -z "$STUDY_ID" ]] && continue  # skip blank lines
    ((COUNT++))

    echo "🔍 [$COUNT] Processing study: $STUDY_ID"

    runs=""

    case "$STUDY_ID" in
        PRJNA*|SRP*)
            runs=$(
                esearch -db sra -query "$STUDY_ID" \
                | efetch -format runinfo \
                | sed 's/\r$//' \
                | awk -F, 'NR>1 {print $1}'
            ) || true
            ;;
        PRJEB*|ERP*|EGAS*)
            runs=$(
                curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$STUDY_ID&result=read_run&fields=run_accession&format=tsv" \
                | sed 's/\r$//' \
                | awk 'NR>1 {print $1}'
            ) || true
            ;;
        GSE*)
            runs=$(
                esearch -db gds -query "$STUDY_ID" \
                | elink -target sra \
                | efetch -format runinfo \
                | sed 's/\r$//' \
                | awk -F, 'NR>1 {print $1}'
            ) || true
            ;;
        *)
            echo "❌ Unsupported study type: $STUDY_ID (skipping)"
            continue
            ;;
    esac

    if [[ -z "$runs" ]]; then
        echo "⚠️  No runs found for $STUDY_ID"
    else
        echo "$runs" >> "$OUTPUT_FILE"
    fi

    echo "✅ Done: $STUDY_ID"
    echo
done < "$STUDY_FILE"

echo "🏁 All studies processed!"
