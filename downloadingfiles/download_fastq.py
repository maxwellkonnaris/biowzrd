#!/usr/bin/env python3
"""
download_fastq.py - Download raw FASTQ files from ENA or SRA.

Usage:
    ./download_fastq.py --accession SRR12345678

Options:
    --accession <str>   Required. SRA or ENA run accession to download.
    --help              Show this help message and exit.

This script:
- First checks ENA for FASTQ files (preferred due to direct FTP access).
- If unavailable, falls back to SRA AWS (using `aws s3 cp`).
- Logs success/failure to a log file (`download.log`).
- Retries up to 3 times on failure.
- Ensures dependencies (`wget`, `aws s3`) are installed.

Example:
    ./download_fastq.py --accession ERR9876543
"""

import os
import sys
import subprocess
import argparse
import requests
import time
from shutil import which

# Directories
WORKDIR = os.getcwd()
FASTQ_DIR = os.path.join(WORKDIR, "fastq_files")
LOG_FILE = os.path.join(WORKDIR, "download.log")
os.makedirs(FASTQ_DIR, exist_ok=True)

# Maximum retries for downloads
MAX_RETRIES = 3

# Check if required dependencies are installed
def check_dependency(command):
    """Check if a system command is available."""
    return which(command) is not None

# Function to log success/failure
def log_message(accession, message):
    """Log success or failure messages."""
    with open(LOG_FILE, "a") as log:
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - {accession} - {message}\n")
    print(message)

# Function to download from ENA
def download_from_ena(accession):
    """Try downloading FASTQ files from ENA."""
    ena_api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp"
    
    try:
        response = requests.get(ena_api_url, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        log_message(accession, f"❌ ENA request failed: {e}")
        return False

    if "fastq_ftp" in response.text:
        fastq_urls = response.text.strip().split("\n")[1].split(";")
        if not fastq_urls or fastq_urls[0] == "":
            log_message(accession, "⚠️ No FASTQ URLs found in ENA.")
            return False

        for url in fastq_urls:
            if url.strip():
                ftp_url = f"ftp://{url}"
                for attempt in range(1, MAX_RETRIES + 1):
                    log_message(accession, f"📥 Downloading from ENA: {ftp_url} (Attempt {attempt}/{MAX_RETRIES})")
                    result = subprocess.run(["wget", "-P", FASTQ_DIR, ftp_url], capture_output=True)
                    if result.returncode == 0:
                        log_message(accession, f"✅ Successfully downloaded from ENA: {ftp_url}")
                        return True
                    time.sleep(2)  # Cooldown between retries

    log_message(accession, "❌ ENA FASTQ download failed.")
    return False

# Function to download from SRA AWS
def download_from_sra(accession):
    """Try downloading FASTQ files from SRA AWS."""
    sra_api_url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={accession}"

    try:
        sra_response = requests.get(sra_api_url, timeout=10)
        sra_response.raise_for_status()
    except requests.RequestException as e:
        log_message(accession, f"❌ SRA AWS request failed: {e}")
        return False

    s3_urls = [line for line in sra_response.text.split() if line.startswith("s3://sra-pub-src")]
    if not s3_urls:
        log_message(accession, "⚠️ No SRA AWS links found.")
        return False

    for s3_url in s3_urls:
        for attempt in range(1, MAX_RETRIES + 1):
            log_message(accession, f"📥 Downloading from SRA AWS: {s3_url} (Attempt {attempt}/{MAX_RETRIES})")
            result = subprocess.run(["aws", "s3", "cp", s3_url, FASTQ_DIR], capture_output=True)
            if result.returncode == 0:
                log_message(accession, f"✅ Successfully downloaded from SRA AWS: {s3_url}")
                return True
            time.sleep(2)  # Cooldown between retries

    log_message(accession, "❌ SRA AWS FASTQ download failed.")
    return False

# Main function
def main():
    """Main execution function."""
    # Argument parser
    parser = argparse.ArgumentParser(description="Download raw FASTQ files from ENA or SRA.")
    parser.add_argument("--accession", type=str, required=True, help="SRA or ENA run accession to download.")
    args = parser.parse_args()

    accession = args.accession.strip()
    
    # Input validation
    if not accession:
        print("❌ ERROR: Accession is required. Use --help for usage.")
        sys.exit(1)

    if not check_dependency("wget"):
        print("❌ ERROR: wget is required but not found. Install it first.")
        sys.exit(1)

    if not check_dependency("aws"):
        print("❌ ERROR: aws CLI is required but not found. Install it first.")
        sys.exit(1)

    log_message(accession, f"🚀 Starting download for {accession}")

    # Try downloading from ENA first
    if not download_from_ena(accession):
        log_message(accession, f"🔄 Falling back to SRA AWS for {accession}")
        if not download_from_sra(accession):
            log_message(accession, f"❌ All download attempts failed for {accession}.")
            sys.exit(1)

    log_message(accession, f"🎉 Download completed for {accession}")
    sys.exit(0)

if __name__ == "__main__":
    main()
