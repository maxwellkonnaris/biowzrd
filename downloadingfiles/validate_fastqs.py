#!/usr/bin/env python3
import gzip
import shutil
import sys
import subprocess
import hashlib
import fcntl
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

# Validate first 64 lines of the FASTQ file
def validate_first_64_lines(fq_path):
    try:
        with gzip.open(fq_path, "rt") as f:
            lines = [next(f).strip() for _ in range(64)]
        if len(lines) < 4 or len(lines) % 4 != 0:
            return False, f"{len(lines)} lines is not a multiple of 4"

        for i in range(0, len(lines), 4):
            if not lines[i].startswith("@"):
                return False, f"Line {i+1} does not start with '@'"
            if not lines[i+2].startswith("+"):
                return False, f"Line {i+3} does not start with '+'"
            if len(lines[i+1]) != len(lines[i+3]):
                return False, f"Seq/qual mismatch at record {i//4+1}"

        return True, "[OK]"
    except Exception as e:
        return False, str(e)

# MD5 checksum calculator
def md5sum(path):
    hash_md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# Attempt to repair a FASTQ file using Biopython
def repair_fastq(fq_path):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    fq_path = Path(fq_path)
    tmp_fastq = fq_path.with_suffix(".repaired.fastq")
    tmp_gz = fq_path.with_suffix(".repaired.fastq.gz")

    try:
        # Count original reads before repair attempt
        with gzip.open(fq_path, "rt") as handle:
            original_reads = sum(1 for _ in FastqGeneralIterator(handle))

        # Attempt to salvage valid reads
        with gzip.open(fq_path, "rt") as in_handle, open(tmp_fastq, "wt") as out_handle:
            repaired_reads = SeqIO.write(SeqIO.parse(in_handle, "fastq"), out_handle, "fastq")

        if repaired_reads == 0:
            tmp_fastq.unlink(missing_ok=True)
            return f"[FAIL] {fq_path}: No valid records found (original had {original_reads} reads)"

        # Compress repaired FASTQ
        compressor = shutil.which("pigz") or shutil.which("gzip")
        compress_cmd = [compressor, "-f", str(tmp_fastq)]
        subprocess.run(compress_cmd, check=True)

        if tmp_gz.exists():
            shutil.move(str(tmp_gz), str(fq_path))
            return (f"[REPAIRED] {fq_path}: replaced with {repaired_reads} salvaged reads "
                    f"from {original_reads} total (corrupted reads were discarded)")
        else:
            return f"[FAIL] {fq_path}: Compression failed after repair"

    except Exception as e:
        tmp_fastq.unlink(missing_ok=True)
        tmp_gz.unlink(missing_ok=True)
        return f"[FAIL] {fq_path}: Repair error: {e}"


# Validate checksum if known, append if not, then check/repair
def check_and_fix_fastq(fq_path, checksums_path, known_checksums, checksums_lock):
    fq_path = Path(fq_path)
    fq_name = fq_path.name

    if fq_name in known_checksums:
        expected_md5 = known_checksums[fq_name]
        actual_md5 = md5sum(fq_path)
        if actual_md5 != expected_md5:
            fq_path.unlink(missing_ok=True)
            return f"[BAD_CHECKSUM] {fq_path}: expected {expected_md5}, got {actual_md5} — file deleted"
    else:
        actual_md5 = md5sum(fq_path)
        with open(checksums_lock, "a+") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            with open(checksums_path, "a") as f:
                f.write(f"{actual_md5}  {fq_name}\n")
            fcntl.flock(lk, fcntl.LOCK_UN)

    valid, message = validate_first_64_lines(fq_path)
    if valid:
        return f"[OK] {fq_path}"
    else:
        repair_result = repair_fastq(fq_path)
        return f"[CORRUPT] {fq_path}: {message} => {repair_result}"

# Load known checksums if present
def load_checksums(checksums_path):
    known = {}
    if checksums_path.exists():
        with open(checksums_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    md5, name = parts
                    known[name] = md5
    return known

# Main driver
def main(base_dir="fastq_data", max_workers=16):
    fq_files = list(Path(base_dir).rglob("*.fastq.gz"))
    if not fq_files:
        print("No FASTQ files found.")
        return

    checksums_path = Path(base_dir) / "checksums.md5"
    checksums_lock = str(checksums_path) + ".lock"
    known_checksums = load_checksums(checksums_path)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                check_and_fix_fastq,
                str(fq),
                checksums_path,
                known_checksums,
                checksums_lock
            ): fq for fq in fq_files
        }
        for fut in as_completed(futures):
            print(fut.result())

if __name__ == "__main__":
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "fastq_data"
    main(base_dir=base_dir, max_workers=16)
