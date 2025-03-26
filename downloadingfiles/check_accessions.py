#!/usr/bin/env python3
import sys
from pathlib import Path
import re

def extract_accessions(filenames):
    accs = set()
    for f in filenames:
        match = re.match(r"(.+?)(_1|_2)?\.fastq\.gz$", f.name)
        if match:
            accs.add(match.group(1))
    return accs

def main():
    run_file = Path("run_accessions.txt")
    if not run_file.exists():
        print("Missing run_accessions.txt")
        sys.exit(1)

    expected = {line.strip() for line in run_file.read_text().splitlines() if line.strip()}
    print(f"📋 Expected: {len(expected)} accessions")

    fastq_dir = Path("fastq_data")
    found = extract_accessions(list(fastq_dir.glob("*.fastq.gz")))

    missing = sorted(expected - found)
    if missing:
        Path("incomplete_run_accessions.txt").write_text("\n".join(missing) + "\n")
        print(f"⚠Missing: {len(missing)} accessions → incomplete_run_accessions.txt")
    else:
        print("All accessions present!")

if __name__ == "__main__":
    main()
