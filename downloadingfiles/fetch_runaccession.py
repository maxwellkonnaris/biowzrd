#!/usr/bin/env python3
import argparse
import csv
import sys
import time
import requests

# Check if Biopython is installed.
try:
    from Bio import Entrez
except ImportError:
    print("Error: Biopython is required but not installed. Please install it (e.g., pip install biopython) and try again.")
    sys.exit(1)

# Attempt to import pysradb for fallback on GEO accessions.
try:
    from pysradb import SRAweb
except ImportError:
    SRAweb = None

def safe_entrez_request(func, *args, max_retries=2, **kwargs):
    for attempt in range(max_retries):
        try:
            handle = func(*args, **kwargs)
            return handle
        except Exception as e:
            if "HTTP Error 429" in str(e):
                wait_time = 2 ** attempt
                print(f"Entrez 429: Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise e
    raise RuntimeError(f"Entrez request failed after {max_retries} retries.")

def safe_requests_get(url, params=None, max_retries=2, timeout=30):
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException as e:
            wait_time = 2 ** attempt
            print(f"Request error: {e}. Retrying in {wait_time}s...")
            time.sleep(wait_time)
    raise RuntimeError(f"Requests failed after {max_retries} retries.")

def get_run_accessions_ncbi(accession, email, api_key):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        handle = safe_entrez_request(
            Entrez.esearch,
            db="sra",
            term=f"{accession}[BioProject]",
            retmax=100000
        )
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            return []

        runs = []
        batch_size = 200
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            handle = safe_entrez_request(
                Entrez.efetch,
                db="sra",
                id=",".join(batch_ids),
                rettype="runinfo",
                retmode="text"
            )
            text = handle.read()
            handle.close()

            if isinstance(text, bytes):
                text = text.decode("utf-8")
            lines = text.strip().split("\n")
            if len(lines) < 2:
                continue

            header = lines[0].split(",")
            try:
                run_index = header.index("Run")
            except ValueError:
                run_index = 0

            for line in lines[1:]:
                cols = line.split(",")
                if len(cols) > run_index:
                    runs.append(cols[run_index].strip())

        return list(set(runs))

    except Exception as e:
        print(f"Error using Bio.Entrez for accession {accession}: {e}")
        return []

def get_run_accessions_geo(accession, email, api_key):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    runs = []

    try:
        handle = safe_entrez_request(Entrez.esearch, db="gds", term=accession, retmax=100000)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if id_list:
            handle = safe_entrez_request(Entrez.elink, dbfrom="gds", db="sra", id=",".join(id_list))
            record = Entrez.read(handle)
            handle.close()

            sra_ids = []
            for rec in record:
                if "LinkSetDb" in rec and rec["LinkSetDb"]:
                    for link in rec["LinkSetDb"][0]["Link"]:
                        sra_ids.append(link["Id"])

            if sra_ids:
                handle = safe_entrez_request(Entrez.efetch, db="sra", id=",".join(sra_ids),
                                             rettype="runinfo", retmode="text")
                text = handle.read()
                handle.close()

                if isinstance(text, bytes):
                    text = text.decode("utf-8")
                lines = text.strip().split("\n")
                if len(lines) >= 2:
                    header = lines[0].split(",")
                    try:
                        run_index = header.index("Run")
                    except ValueError:
                        run_index = 0
                    for line in lines[1:]:
                        cols = line.split(",")
                        if len(cols) > run_index:
                            runs.append(cols[run_index].strip())
    except Exception as e:
        print(f"Entrez failed for GEO accession {accession}: {e}")

    if not runs and SRAweb is not None:
        db = SRAweb()
        for attempt in range(2):
            try:
                print(f"Falling back to pysradb for GEO accession {accession} (attempt {attempt + 1})")
                df = db.sra_metadata(geo=accession, detailed=True)
                if not df.empty and "run_accession" in df.columns:
                    runs.extend(df["run_accession"].dropna().unique().tolist())
                break
            except Exception as e:
                wait_time = 2 ** attempt
                print(f"pysradb error: {e}. Retrying in {wait_time}s...")
                time.sleep(wait_time)

    return list(set(runs))

def fetch_runs_ena(accession):
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession",
        "format": "tsv"
    }
    try:
        r = safe_requests_get(url, params=params)
        lines = r.text.strip().split("\n")
        if len(lines) > 1:
            return list(set([line.strip() for line in lines[1:] if line.strip()]))
        else:
            return []
    except Exception as e:
        print(f"ENA fetch failed for {accession}: {e}")
        return []

def process_accession(accession, email, api_key):
    accession = accession.strip()
    if accession.startswith("PRJNA") or accession.startswith("SRP"):
        runs_ncbi = get_run_accessions_ncbi(accession, email, api_key)
        runs_ena = fetch_runs_ena(accession)
        return list(set(runs_ncbi + runs_ena))
    elif accession.startswith("GSE"):
        return get_run_accessions_geo(accession, email, api_key)
    elif accession.startswith("PRJEB") or accession.startswith("ERP") or accession.startswith("EGAS"):
        return fetch_runs_ena(accession)
    else:
        print(f"Unknown accession prefix for {accession}. Skipping.")
        return []

def read_accessions(input_file):
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
    parser.add_argument("input_file", help="Input file with one project accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="CSV output file name.")
    parser.add_argument("--runs-only", default="run_accessions_only.txt", help="Text file with only run accessions.")
    parser.add_argument("--fail-log", default="failed_accessions.log", help="File to log failed accessions.")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    args = parser.parse_args()

    email = args.email if args.email else input("Enter your email (required for NCBI queries): ").strip()
    api_key = args.api_key if args.api_key else input("Enter your NCBI API key (press enter if none): ").strip()

    if not email:
        print("Error: Email is required for NCBI queries.")
        sys.exit(1)

    accessions = read_accessions(args.input_file)
    results = []
    failed = []
    runs_only = []
    seen_runs = set()

    for accession in accessions:
        print(f"\n🔍 Processing: {accession}")
        try:
            run_list = process_accession(accession, email, api_key)
            if run_list:
                for run in run_list:
                    if run not in seen_runs:
                        results.append((accession, run))
                        runs_only.append(run)
                        seen_runs.add(run)
            else:
                failed.append((accession, "No run accessions found"))
        except Exception as e:
            print(f"❌ Failed to process {accession}: {e}")
            failed.append((accession, str(e)))
        time.sleep(0.5)

    try:
        with open(args.output_file, "w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Project Accession", "Run Accession"])
            writer.writerows(results)
        print(f"\n✅ Results saved to {args.output_file}")
    except Exception as e:
        print(f"Error writing CSV: {e}")

    try:
        with open(args.runs_only, "w") as txtfile:
            for run in runs_only:
                txtfile.write(run + "\n")
        print(f"✅ Run accessions only saved to {args.runs_only}")
    except Exception as e:
        print(f"Error writing run accessions only file: {e}")

    if failed:
        try:
            with open(args.fail_log, "w") as flog:
                for acc, reason in failed:
                    flog.write(f"{acc}\t{reason}\n")
            print(f"⚠️ Failed accessions logged to {args.fail_log}")
        except Exception as e:
            print(f"Error writing failure log: {e}")

if __name__ == '__main__':
    main()
