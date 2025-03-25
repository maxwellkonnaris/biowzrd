#!/usr/bin/env python3
import argparse
import csv
import re
import sys
import time
import requests

# Check if Biopython is installed.
try:
    from Bio import Entrez
except ImportError:
    print("Error: Biopython is required but not installed. Please install it (e.g., pip install biopython) and try again.")
    sys.exit(1)

def get_run_accessions_ncbi(accession, email, api_key):
    """
    For accessions starting with PRJNA* or SRP*, use Bio.Entrez to search SRA,
    fetch run information, and return a list of run accessions.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    try:
        # Search SRA for the accession
        handle = Entrez.esearch(db="sra", term=accession)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            return []
        # Fetch run info using the list of IDs.
        try:
            handle = Entrez.efetch(db="sra", id=",".join(id_list), rettype="runinfo", retmode="text")
            text = handle.read()
            handle.close()
        except Exception as e:
            if "HTTP Error 429" in str(e):
                print(f"HTTP 429 received for {accession}, waiting before retrying...")
                time.sleep(5)
                handle = Entrez.efetch(db="sra", id=",".join(id_list), rettype="runinfo", retmode="text")
                text = handle.read()
                handle.close()
            else:
                raise e
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        lines = text.strip().split("\n")
        if len(lines) < 2:
            return []
        header = lines[0].split(",")
        try:
            run_index = header.index("Run")
        except ValueError:
            run_index = 0  # fallback if "Run" is not found
        runs = []
        for line in lines[1:]:
            cols = line.split(",")
            if len(cols) > run_index:
                runs.append(cols[run_index].strip())
        return runs
    except Exception as e:
        print(f"Error using Bio.Entrez for accession {accession}: {e}")
        return []

def get_run_accessions_geo(accession, email, api_key):
    """
    For GEO accessions (GSE*), search the gds database,
    then use elink to connect to SRA and fetch run info.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    try:
        # Search gds for the GEO accession
        handle = Entrez.esearch(db="gds", term=accession)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            return []
        # Use elink to link from gds to sra
        handle = Entrez.elink(dbfrom="gds", db="sra", id=",".join(id_list))
        record = Entrez.read(handle)
        handle.close()
        sra_ids = []
        for rec in record:
            if "LinkSetDb" in rec and rec["LinkSetDb"]:
                for link in rec["LinkSetDb"][0]["Link"]:
                    sra_ids.append(link["Id"])
        if not sra_ids:
            return []
        # Fetch run info using the linked SRA IDs.
        try:
            handle = Entrez.efetch(db="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text")
            text = handle.read()
            handle.close()
        except Exception as e:
            if "HTTP Error 429" in str(e):
                print(f"HTTP 429 received for GEO accession {accession}, waiting before retrying...")
                time.sleep(5)
                handle = Entrez.efetch(db="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text")
                text = handle.read()
                handle.close()
            else:
                raise e
        if isinstance(text, bytes):
            text = text.decode("utf-8")
        lines = text.strip().split("\n")
        if len(lines) < 2:
            return []
        header = lines[0].split(",")
        try:
            run_index = header.index("Run")
        except ValueError:
            run_index = 0
        runs = []
        for line in lines[1:]:
            cols = line.split(",")
            if len(cols) > run_index:
                runs.append(cols[run_index].strip())
        return runs
    except Exception as e:
        print(f"Error using Bio.Entrez for GEO accession {accession}: {e}")
        return []

def fetch_runs_ena(accession):
    """
    For accessions starting with PRJEB*, ERP*, or EGAS*, use the ENA API to fetch run accessions.
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
    Select the appropriate method for fetching run accessions based on the accession prefix.
    """
    accession = accession.strip()
    runs = []
    if accession.startswith("PRJNA") or accession.startswith("SRP"):
        runs = get_run_accessions_ncbi(accession, email, api_key)
    elif accession.startswith("GSE"):
        runs = get_run_accessions_geo(accession, email, api_key)
    elif accession.startswith("PRJEB") or accession.startswith("ERP") or accession.startswith("EGAS"):
        runs = fetch_runs_ena(accession)
    else:
        print(f"Unknown accession prefix for {accession}. Skipping.")
    return runs

def read_accessions(input_file):
    """
    Read the input file and return a list of project accessions.
    The file should contain one accession per line (ignores blank lines and comments starting with '#').
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
        epilog="Example usage: python download_runs.py projects.txt --email your_email@example.com --api-key YOUR_API_KEY"
    )
    parser.add_argument("input_file", help="Input file (txt/tsv/csv) with one project accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="Output CSV file name (default: run_accessions.csv).")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    args = parser.parse_args()

    # Get email and API key (prompt if not provided via flags)
    email = args.email if args.email else input("Enter your email (required for NCBI queries): ").strip()
    api_key = args.api_key if args.api_key else input("Enter your NCBI API key (press enter if none): ").strip()

    accessions = read_accessions(args.input_file)
    results = []  # List of tuples: (Project Accession, Run Accession)

    for accession in accessions:
        print(f"Processing accession: {accession}")
        run_list = process_accession(accession, email, api_key)
        if run_list:
            for run in run_list:
                results.append((accession, run))
        else:
            results.append((accession, "No run accessions found"))
        # Pause briefly to avoid overloading servers
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
