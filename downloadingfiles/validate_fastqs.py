#!/usr/bin/env python3
import gzip
import shutil
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

def validate_first_64_lines(fq_path):
    try:
        with gzip.open(fq_path, "rt") as f:
            lines = [next(f).strip() for _ in range(64)]
        if len(lines) < 4 or len(lines) % 4 != 0:
            return False, f"{len(lines)} lines is not valid FASTQ"
        for i in range(0, len(lines), 4):
            if not lines[i].startswith("@"):
                return False, f"Line {i+1} does not start with '@'"
            if not lines[i+2].startswith("+"):
                return False, f"Line {i+3} does not start with '+'"
            if len(lines[i+1]) != len(lines[i+3]):
                return False, f"Seq/qual length mismatch at record {i//4+1}"
        return True, "[OK]"
    except Exception as e:
        return False, str(e)

def repair_fastq(fq_path):
    fq_path = Path(fq_path)
    tmp_fastq = fq_path.with_suffix(".repaired.fastq")
    tmp_gz = fq_path.with_suffix(".repaired.fastq.gz")

    try:
        # Step 1: Parse and write uncompressed FASTQ
        with gzip.open(fq_path, "rt") as in_handle, open(tmp_fastq, "wt") as out_handle:
            count = SeqIO.write(SeqIO.parse(in_handle, "fastq"), out_handle, "fastq")

        if count == 0:
            tmp_fastq.unlink(missing_ok=True)
            return f"[FAIL] {fq_path}: No records found in repair"

        # Step 2: Compress with pigz or gzip
        compressor = shutil.which("pigz") or shutil.which("gzip")
        compress_cmd = [compressor, "-f", str(tmp_fastq)]
        subprocess.run(compress_cmd, check=True)

        # Step 3: Rename repaired .gz to original
        if tmp_gz.exists():
            shutil.move(str(tmp_gz), str(fq_path))
            return f"[REPAIRED] {fq_path} ({count} reads)"
        else:
            return f"[FAIL] {fq_path}: Compression failed"

    except Exception as e:
        tmp_fastq.unlink(missing_ok=True)
        tmp_gz.unlink(missing_ok=True)
        return f"[FAIL] {fq_path}: Repair error: {e}"

def check_and_fix_fastq(fq_path):
    fq_path = Path(fq_path)
    valid, message = validate_first_64_lines(fq_path)
    if valid:
        return f"[OK] {fq_path}"
    else:
        repair_result = repair_fastq(fq_path)
        return f"[CORRUPT] {fq_path}: {message} => {repair_result}"

def main(base_dir="fastq_data", max_workers=16):
    fq_files = list(Path(base_dir).rglob("*.fastq.gz"))
    if not fq_files:
        print("No FASTQ files found.")
        return

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(check_and_fix_fastq, str(fq)): fq for fq in fq_files}
        for fut in as_completed(futures):
            print(fut.result())

if __name__ == "__main__":
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "fastq_data"
    main(base_dir=base_dir, max_workers=16)
