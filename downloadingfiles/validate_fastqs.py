#!/usr/bin/env python3
import gzip
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def validate_first_64_lines(fq_path):
    try:
        with gzip.open(fq_path, "rt") as f:
            lines = []
            for i, line in enumerate(f):
                lines.append(line.strip())
                if i == 63:
                    break

        if len(lines) < 4:
            return f"[FAIL] {fq_path}: Too few lines ({len(lines)})"

        if len(lines) % 4 != 0:
            return f"[FAIL] {fq_path}: Not multiple of 4 lines ({len(lines)})"

        for i in range(0, len(lines), 4):
            if not lines[i].startswith("@"):
                return f"[FAIL] {fq_path}: Line {i+1} does not start with '@'"
            if not lines[i+2].startswith("+"):
                return f"[FAIL] {fq_path}: Line {i+3} does not start with '+'"
            if len(lines[i+1]) != len(lines[i+3]):
                return f"[FAIL] {fq_path}: Seq/qual length mismatch (line {i+2})"

        return f"[OK] {fq_path}"

    except Exception as e:
        return f"[FAIL] {fq_path}: {e}"

def main(base_dir="fastq_data", max_workers=8):
    fq_files = list(Path(base_dir).rglob("*.fastq.gz"))
    if not fq_files:
        print("No FASTQ files found.")
        return

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(validate_first_64_lines, str(fq)): fq for fq in fq_files}
        for fut in as_completed(futures):
            result = fut.result()
            print(result)

if __name__ == "__main__":
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "fastq_data"
    main(base_dir=base_dir, max_workers=16)
