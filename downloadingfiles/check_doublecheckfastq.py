#!/usr/bin/env python3
import os
import subprocess
import sys
import glob
from pathlib import Path
import shutil

REQUIRED_TOOLS = ["fasterq-dump", "pigz"]

def check_tools_exist(tools):
    """Ensure required tools are available in PATH."""
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        print(f"ERROR: Missing required tools in PATH: {', '.join(missing)}")
        sys.exit(1)

def run(cmd, error_msg):
    """Run a command with error handling."""
    try:
        print(f"🔹 Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"{error_msg}: {e}")
        sys.exit(1)

def main():
    check_tools_exist(REQUIRED_TOOLS)

    array_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if array_id is None:
        print("SLURM_ARRAY_TASK_ID not set. This must be run as a Slurm array job.")
        sys.exit(1)

    fastq_dir = Path("./fastq_data")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    sra_list_path = Path("sra_list.txt")
    if not sra_list_path.exists():
        print("sra_list.txt not found in current directory.")
        sys.exit(1)

    with sra_list_path.open() as f:
        sra_list = [line.strip() for line in f if line.strip()]

    try:
        sra_file = sra_list[int(array_id) - 1]
    except (IndexError, ValueError):
        print(f"Invalid SLURM_ARRAY_TASK_ID: {array_id}")
        sys.exit(1)

    sra_path = fastq_dir / sra_file
    if not sra_path.exists():
        print(f"Expected SRA file not found: {sra_path}")
        sys.exit(1)

    basename = sra_path.stem
    print(f"[Task {array_id}] Starting on {sra_path.name}")

    # Run fasterq-dump
    run([
        "fasterq-dump",
        str(sra_path),
        "--outdir", str(fastq_dir),
        "--threads", "4",
        "--mem", "8G",
        "--split-3"
    ], f"[Task {array_id}] fasterq-dump failed")

    # Compress FASTQ files using pigz
    fastq_files = glob.glob(str(fastq_dir / f"{basename}*.fastq"))
    if not fastq_files:
        print(f"No FASTQ files found to compress for {basename}")
        sys.exit(1)

    pigz_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    compress_success = True

    for fq in fastq_files:
        try:
            run(["pigz", "-p", str(pigz_threads), fq], f"[Task {array_id}] pigz failed on {fq}")
        except Exception as e:
            print(f"Compression error: {e}")
            compress_success = False

    if compress_success:
        print(f"[Task {array_id}] Compression complete. Deleting {sra_path.name}")
        try:
            sra_path.unlink()
        except Exception as e:
            print(f"Could not delete {sra_path}: {e}")
    else:
        print(f"[Task {array_id}] Compression failed. Keeping {sra_path.name} for inspection.")

if __name__ == "__main__":
    main()
