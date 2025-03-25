#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import time
import shutil
import sys
import requests

def run_ncbi_pipeline(accession, email, api_key, database, pipeline):
    """
    Run a pipeline using NCBI's esearch/efetch (and optionally elink) tools.
    Builds a command string that queries the given database with esearch (using email and API key),
    then pipes the output through the provided pipeline and returns a list of run accessions.
    """
    # Check for required command-line tools.
    tools_needed = ["esearch", "efetch"]
    if "elink" in pipeline:
        tools_needed.append("elink")
    for tool in tools_needed:
        if shutil.which(tool) is None:
            print(f"Warning: {tool} not found in PATH. Please install NCBI EDirect tools if you wish to use the NCBI query mode.")
            return []

    # Build the base command.
    cmd = f"esearch -db {database} -query {accession} -email {email}"
    if api_key:
        cmd += f" -api_key {api_key}"
    # Append the additional pipeline commands (efetch, sed, awk, etc.)
    cmd += pipeline

    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
        # Split the output into lines and return non-empty ones.
        lines = output.strip().split("\n")
        runs = [line.strip() for line in lines if line.strip()]
        return runs
    except subprocess.CalledProcessError as e:
        print(f"Error running command for accession {accession}: {e}")
        return []

def fetch_runs_ena(accession):
    """
    For accessions that match PRJEB*, ERP* or EGAS* patterns,
    use the ENA API to fetch run accessions.
    """
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession",
        "format": "tsv"
    }
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        # Skip header if present.
        if len(lines) > 1:
            runs = [line.strip() for line in lines[1:] if line.strip()]
        else:
            runs = []
        return runs
    except Exception as e:
        print(f"Error fetching ENA data for {accession}: {e}")
        return []

def process_accession(accession, email, api_key):
    """
    Choose the appropriate method for fetching run accessions based on the accession pattern.
    """
    accession = accession.strip()
    runs = []
    if re.match(r'^(PRJNA|SRP)', accession):
        # Use NCBI EDirect commands for SRA accessions.
        pipeline = " | efetch -format runinfo | sed 's/\\r$//' | awk -F, 'NR>1 {print $1}'"
        runs = run_ncbi_pipeline(accession, email, api_key, database="sra", pipeline=pipeline)
    elif re.match(r'^(PRJEB|ERP|EGAS)', accession):
        # Use the ENA API.
        runs = fetch_runs_ena(accession)
    elif accession.startswith("GSE"):
        # For GEO accessions, use esearch on the gds database, then link to SRA.
        pipeline = " | elink -target sra | efetch -format runinfo | sed 's/\\r$//' | awk -F, 'NR>1 {print $1}'"
        runs = run_ncbi_pipeline(accession, email, api_key, database="gds", pipeline=pipeline)
    else:
        print(f"Unknown accession prefix for {accession}. Skipping.")
    return runs

def read_accessions(input_file):
    """
    Read the input file and return a list of project accessions.
    Assumes one accession per line (ignores blank lines and lines starting with #).
    """
    accessions = []
    try:
        with open(input_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    accessions.append(line)
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        sys.exit(1)
    return accessions

def main():
    parser = argparse.ArgumentParser(
        description="Download run accessions for a list of project accessions sequentially.",
        epilog="Example usage: python download_runs.py projects.txt --email user@example.com --api-key YOUR_API_KEY"
    )
    parser.add_argument("input_file", help="Input file (txt/tsv/csv) with one project accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="Output CSV file name (default: run_accessions.csv).")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    args = parser.parse_args()

    # Check for email; if not provided, prompt the user.
    email = args.email if args.email else input("Enter your email (required for NCBI queries): ").strip()
    api_key = args.api_key if args.api_key else input("Enter your NCBI API key (press enter if none): ").strip()

    accessions = read_accessions(args.input_file)
    results = []  # Will store tuples: (Project Accession, Run Accession)

    for accession in accessions:
        print(f"Processing accession: {accession}")
        run_list = process_accession(accession, email, api_key)
        if run_list:
            for run in run_list:
                results.append((accession, run))
        else:
            results.append((accession, "No run accessions found"))
        # Sequential processing with a brief pause to avoid overloading the servers.
        time.sleep(0.5)

    # Write the results sequentially to a CSV file.
    try:
        with open(args.output_file, mode="w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Project Accession", "Run Accession"])
            writer.writerows(results)
        print(f"Results saved to {args.output_file}")
    except Exception as e:
        print(f"Error writing to {args.output_file}: {e}")

if __name__ == '__main__':
    main()
