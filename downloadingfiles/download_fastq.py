#!/usr/bin/env python3
"""
download_fastq.py

Download raw FASTQ files in parallel from paths specified in a tab-separated file.
Each line must contain two columns: 'accession' and 'download_path'.

Example input (download_links.txt):

    accession       download_path
    SRR32578126     s3://sra-pub-src-14/SRR32578126/file_R1.fastq.gz
    SRR32578126     s3://sra-pub-src-14/SRR32578126/file_R2.fastq.gz
    SRR99999999     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR999/...
    ...

Usage:
    python download_fastq.py --input download_links.txt [--output raw_fastq_files] [--verbose]

Options:
    --input <str>       Path to the tab-separated input file (required).
    --output <str>      Directory to store downloaded FASTQ files (default: ./raw_fastq_files).
    --verbose           Print log messages to console in addition to writing them to download.log.
    --help              Show this help message and exit.

Details:
- Script attempts to detect s3:// vs. any other protocol (ftp/http/https).
- Retries up to MAX_RETRIES times on failure.
- Uses 90% of detected CPU cores for parallel downloads.
- Writes SUCCESS/FAILURE to download.log, and to console if --verbose is used.
"""

import os
import sys
import time
import argparse
import subprocess
from shutil import which
from multiprocessing import Pool, cpu_count

# Globals
MAX_RETRIES = 3
LOG_FILE = "download.log"

def log_message(msg, log_file=LOG_FILE, verbose=False):
    """
    Log a message to file and also to console if verbose is True.
    """
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    full_msg = f"{timestamp} - {msg}"
    # Always log to file
    with open(log_file, "a") as lf:
        lf.write(full_msg + "\n")
    # Print to console if verbose
    if verbose:
        print(full_msg)

def check_dependency(command):
    """
    Check if a system command is available. Return True if found, otherwise False.
    """
    return which(command) is not None

def download_single(accession, download_path, output_dir, verbose=False):
    """
    Download a single file (accession + download_path).
    - If download_path starts with 's3://', use `aws s3 cp`.
    - Otherwise, use `wget`.
    Logs SUCCESS or FAILURE for each line (accession).
    Retries up to MAX_RETRIES times.
    """
    # Ensure output_dir exists
    os.makedirs(output_dir, exist_ok=True)

    # Decide how to download based on prefix
    if download_path.startswith("s3://"):
        # Make sure 'aws' is installed
        if not check_dependency("aws"):
            log_message(f"{accession} - FAILURE: 'aws' CLI not found.", verbose=verbose)
            return
        command = ["aws", "s3", "cp", download_path, output_dir]
    else:
        # Make sure 'wget' is installed
        if not check_dependency("wget"):
            log_message(f"{accession} - FAILURE: 'wget' not found.", verbose=verbose)
            return
        command = ["wget", "-P", output_dir, download_path]

    success = False
    for attempt in range(1, MAX_RETRIES + 1):
        log_message(f"{accession} - Attempt {attempt}/{MAX_RETRIES} - Downloading {download_path}", verbose=verbose)
        result = subprocess.run(command, capture_output=True)
        if result.returncode == 0:
            # Download succeeded
            log_message(f"{accession} - SUCCESS: Downloaded {download_path}", verbose=verbose)
            success = True
            break
        else:
            # Download failed
            log_message(
                f"{accession} - WARNING: Failed to download {download_path}\n"
                f"stderr: {result.stderr.decode('utf-8', errors='replace')}",
                verbose=verbose
            )
            time.sleep(2)  # short cooldown before next retry

    if not success:
        log_message(f"{accession} - FAILURE: Could not download {download_path} after {MAX_RETRIES} attempts", 
                    verbose=verbose)

def main():
    parser = argparse.ArgumentParser(
        description="Download FASTQ files from a tab-separated file of 'accession' and 'download_path'."
    )
    parser.add_argument("--input", required=True, help="Path to the tab-separated input file.")
    parser.add_argument("--output", default="raw_fastq_files", help="Output directory for downloaded files.")
    parser.add_argument("--verbose", action="store_true", help="Print log messages to console as well.")
    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output
    verbose = args.verbose

    # Determine how many processes to run
    num_cores = cpu_count()
    num_processes = max(1, int(num_cores * 0.9))

    # Read the input file lines
    tasks = []
    with open(input_file, "r") as f:
        # Skip header if present, or handle it safely
        header = f.readline().strip().split()
        # A naive check if it's the right header; you can adapt as needed.
        if "accession" in header and "download_path" in header:
            # We already consumed the header line; move on
            pass
        else:
            # If no header recognized, then treat that line as data
            # Rewind or parse differently
            f.seek(0)

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                # Not a well‐formed line
                continue
            accession, download_path = parts[0], parts[1]
            tasks.append((accession, download_path, output_dir, verbose))

    if not tasks:
        print("No valid lines found in input file. Exiting.")
        sys.exit(1)

    log_message(f"Starting downloads using {num_processes} parallel processes.", verbose=verbose)

    # Use a multiprocessing Pool to parallelize
    with Pool(processes=num_processes) as pool:
        pool.starmap(download_single, tasks)

    log_message("All download tasks completed (SUCCESS or FAILURE logged).", verbose=verbose)

if __name__ == "__main__":
    main()
