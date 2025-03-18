import os
import sys
import subprocess
import requests

# Set up directories
WORKDIR = os.getcwd()
FASTQ_DIR = os.path.join(WORKDIR, "fastq_files")
os.makedirs(FASTQ_DIR, exist_ok=True)

# Read accession from SLURM input
if len(sys.argv) != 2:
    print("Usage: python download_fastq.py <ACCESSION>")
    sys.exit(1)

accession = sys.argv[1]
print(f"\nProcessing {accession}...")

def download_from_ena(accession):
    ena_api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp"
    response = requests.get(ena_api_url)

    if response.status_code == 200 and "fastq_ftp" in response.text:
        fastq_urls = response.text.split("\n")[1].split(";")
        for url in fastq_urls:
            if url.strip():
                print(f"Downloading from ENA: {url}")
                subprocess.run(["wget", "-P", FASTQ_DIR, f"ftp://{url}"])
        return True
    return False

def download_from_sra(accession):
    sra_api_url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run={accession}"
    sra_response = requests.get(sra_api_url)

    if "s3://sra-pub-src" in sra_response.text:
        s3_urls = [line for line in sra_response.text.split() if line.startswith("s3://sra-pub-src")]
        for s3_url in s3_urls:
            print(f"Downloading from SRA AWS: {s3_url}")
            subprocess.run(["aws", "s3", "cp", s3_url, FASTQ_DIR])
        return True
    return False

if not download_from_ena(accession):
    print(f"FASTQ not found in ENA. Checking SRA AWS for {accession}...")
    if not download_from_sra(accession):
        print(f"FASTQ files not found for {accession}. Check manually.")

print("Download complete!")
