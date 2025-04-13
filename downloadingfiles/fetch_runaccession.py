#!/usr/bin/env python3
import argparse
import csv
import sys
import time
import requests
import re
import xml.etree.ElementTree as ET
import json
import signal
import os
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

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

# Attempt to import tqdm for progress bar.
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
    print("Warning: tqdm not installed. Progress bar will be disabled in single-threaded mode. Install with 'pip install tqdm'.")

# Global lock for thread-safe operations
print_lock = Lock()
checkpoint_lock = Lock()

# Global checkpoint data
checkpoint_data = {
    "completed": [],
    "results": [],  # List of (accession, run) tuples
    "failed": [],   # List of (accession, reason) tuples
    "comparison": [],  # List of {accession, expected, found, status}
    "runs_only": [],   # List of run accessions
}

def save_checkpoint():
    """Save checkpoint data to checkpoint.json."""
    global checkpoint_data
    with checkpoint_lock:
        try:
            with open("checkpoint.json", "w") as f:
                json.dump(checkpoint_data, f, indent=2)
        except Exception as e:
            print(f"Error saving checkpoint: {e}")

def load_checkpoint():
    """Load checkpoint data from checkpoint.json if it exists."""
    global checkpoint_data
    if os.path.exists("checkpoint.json"):
        try:
            with open("checkpoint.json", "r") as f:
                checkpoint_data.update(json.load(f))
            print(f"Loaded checkpoint: {len(checkpoint_data['completed'])} accessions already processed.")
        except Exception as e:
            print(f"Error loading checkpoint: {e}")

def signal_handler(sig, frame):
    """Handle interrupt (Ctrl+C) by saving checkpoint and exiting."""
    print("\nInterrupt received, saving checkpoint...")
    save_checkpoint()
    print("Checkpoint saved to checkpoint.json. Resume by rerunning the script.")
    sys.exit(0)

def safe_entrez_request(func, *args, max_retries=2, verbose=False, **kwargs):
    """Helper to gracefully retry Entrez calls on HTTP 429 or transient failures."""
    for attempt in range(max_retries):
        try:
            handle = func(*args, **kwargs)
            return handle
        except Exception as e:
            if "HTTP Error 429" in str(e):
                wait_time = 2 ** attempt
                if verbose:
                    with print_lock:
                        print(f"Entrez 429: Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                raise e
    raise RuntimeError(f"Entrez request failed after {max_retries} retries.")

def safe_requests_get(url, params=None, max_retries=2, timeout=30, verbose=False):
    """Helper to retry requests.get on transient failures."""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.exceptions.HTTPError as e:
            wait_time = 2 ** attempt
            if e.response.status_code == 400:
                if verbose:
                    with print_lock:
                        print(f"Bad request (400) for {url}: {e}. Retrying in {wait_time}s...")
            else:
                if verbose:
                    with print_lock:
                        print(f"Request error: {e}. Retrying in {wait_time}s...")
            time.sleep(wait_time)
        except requests.exceptions.RequestException as e:
            wait_time = 2 ** attempt
            if verbose:
                with print_lock:
                    print(f"Request error: {e}. Retrying in {wait_time}s...")
            time.sleep(wait_time)
    raise RuntimeError(f"Requests failed after {max_retries} retries.")

def fetch_runs_ena(accession, verbose=False):
    """Query ENA Filereport for run_accession."""
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession",
        "format": "tsv"
    }
    try:
        r = safe_requests_get(url, params=params, verbose=verbose)
        lines = r.text.strip().split("\n")
        if len(lines) > 1:
            return list(set(line.strip() for line in lines[1:] if line.strip()))
        return []
    except RuntimeError as e:
        if "400 Client Error" in str(e):
            if verbose:
                with print_lock:
                    print(f"ENA fetch failed for {accession}: Invalid accession or no runs available.")
        else:
            if verbose:
                with print_lock:
                    print(f"ENA fetch failed for {accession}: {e}")
        return []

def get_run_accessions_ncbi(accession, email, api_key, verbose=False):
    """Fetch run accessions from NCBI SRA for BioProject (PRJNA/SRP), including SRS and XML RUN_SET."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        # Step 1: Search SRA for BioProject
        if verbose:
            with print_lock:
                print(f" - Searching SRA for BioProject {accession}...")
        handle = safe_entrez_request(
            Entrez.esearch,
            db="sra",
            term=f"{accession}[BioProject]",
            retmax=100000,
            verbose=verbose
        )
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            if verbose:
                with print_lock:
                    print(f" - No SRA experiments found for {accession}.")
            return [], 0

        # Initial estimate of expected runs
        expected_runs = len(id_list)  # Proxy: one run per SRX
        runs = []
        batch_size = 200
        # Step 2: Fetch runinfo for initial IDs
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            handle = safe_entrez_request(
                Entrez.efetch,
                db="sra",
                id=",".join(batch_ids),
                rettype="runinfo",
                retmode="text",
                verbose=verbose
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

            # Collect runs from runinfo
            for line in lines[1:]:
                cols = line.split(",")
                if len(cols) > run_index and cols[run_index].strip():
                    runs.append(cols[run_index].strip())

        # Step 3: Check XML for SRS and runs using ElementTree
        if verbose:
            with print_lock:
                print(f" - Checking XML for SRS and runs in {accession}...")
        srs_ids = []
        xml_runs = []
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            handle = safe_entrez_request(
                Entrez.efetch,
                db="sra",
                id=",".join(batch_ids),
                rettype="full",
                retmode="xml",
                verbose=verbose
            )
            xml_text = handle.read()
            handle.close()

            if isinstance(xml_text, bytes):
                xml_text = xml_text.decode("utf-8")

            try:
                root = ET.fromstring(xml_text)
                # Parse XML for SRS and runs
                for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
                    # Extract SRS
                    sample = exp_pkg.find(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
                    if sample is not None and sample.text and sample.text.startswith("SRS"):
                        srs_ids.append(sample.text)
                    # Extract runs from RUN_SET
                    run_set = exp_pkg.find(".//RUN_SET")
                    if run_set is not None:
                        for run in run_set.findall("RUN"):
                            run_acc = run.get("accession")
                            if run_acc and run_acc.startswith(("SRR", "ERR", "DRR")):
                                xml_runs.append(run_acc)
            except ET.ParseError as e:
                with print_lock:
                    print(f"XML parsing failed for {accession}: {e}")
                continue

        if xml_runs:
            if verbose:
                with print_lock:
                    print(f" - Found runs in XML RUN_SET for {accession}.")
            runs.extend(xml_runs)
            # Update expected runs based on XML RUN_SET
            expected_runs = max(expected_runs, len(set(xml_runs)))

        # Step 4: Query SRA for runs linked to SRS
        if srs_ids:
            if verbose:
                with print_lock:
                    print(f" - Found SRS for {accession}. Fetching runs...")
            for srs in srs_ids:
                handle = safe_entrez_request(
                    Entrez.esearch,
                    db="sra",
                    term=f"{srs}[Sample]",
                    retmax=100000,
                    verbose=verbose
                )
                srs_record = Entrez.read(handle)
                handle.close()
                srs_id_list = srs_record.get("IdList", [])
                if srs_id_list:
                    for j in range(0, len(srs_id_list), batch_size):
                        srs_batch_ids = srs_id_list[j:j + batch_size]
                        handle = safe_entrez_request(
                            Entrez.efetch,
                            db="sra",
                            id=",".join(srs_batch_ids),
                            rettype="runinfo",
                            retmode="text",
                            verbose=verbose
                        )
                        srs_text = handle.read()
                        handle.close()

                        if isinstance(srs_text, bytes):
                            srs_text = srs_text.decode("utf-8")
                        srs_lines = srs_text.strip().split("\n")
                        if len(srs_lines) < 2:
                            continue

                        srs_header = srs_lines[0].split(",")
                        try:
                            srs_run_index = srs_header.index("Run")
                        except ValueError:
                            srs_run_index = 0

                        for line in srs_lines[1:]:
                            cols = line.split(",")
                            if len(cols) > srs_run_index and cols[srs_run_index].strip():
                                runs.append(cols[srs_run_index].strip())

        # Final expected runs: max of SRX count, runinfo runs, and XML runs
        expected_runs = max(expected_runs, len(set(runs)))
        return list(set(runs)), expected_runs

    except Exception as e:
        with print_lock:
            print(f"Error processing NCBI accession {accession}: {e}")
        return [], 0

def get_run_accessions_geo(accession, email, api_key, verbose=False):
    """Resolve GSE to SRR."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    runs = []
    try:
        handle = safe_entrez_request(Entrez.esearch, db="gds", term=accession, retmax=100000, verbose=verbose)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if id_list:
            handle = safe_entrez_request(Entrez.elink, dbfrom="gds", db="sra", id=",".join(id_list), verbose=verbose)
            record = Entrez.read(handle)
            handle.close()

            sra_ids = []
            for rec in record:
                if "LinkSetDb" in rec and rec["LinkSetDb"]:
                    for link in rec["LinkSetDb"][0]["Link"]:
                        sra_ids.append(link["Id"])

            if sra_ids:
                handle = safe_entrez_request(Entrez.efetch, db="sra", id=",".join(sra_ids),
                                             rettype="runinfo", retmode="text", verbose=verbose)
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
        with print_lock:
            print(f"Entrez failed for GEO accession {accession}: {e}")

    if not runs and SRAweb is not None:
        db = SRAweb()
        for attempt in range(2):
            try:
                if verbose:
                    with print_lock:
                        print(f"Falling back to pysradb for GEO accession {accession} (attempt {attempt + 1})")
                df = db.sra_metadata(geo=accession, detailed=True)
                if not df.empty and "run_accession" in df.columns:
                    runs.extend(df["run_accession"].dropna().unique().tolist())
                break
            except Exception as e:
                wait_time = 2 ** attempt
                if verbose:
                    with print_lock:
                        print(f"pysradb error: {e}. Retrying in {wait_time}s...")
                time.sleep(wait_time)

    return list(set(runs)), 0  # No expected runs estimate for GEO

def fallback_entrez_accession(accession, email, api_key, verbose=False):
    """Fallback: Search for runs by accession[ACCN], including SRS and XML RUN_SET."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        handle = safe_entrez_request(
            Entrez.esearch,
            db="sra",
            term=f"{accession}[ACCN]",
            retmax=100000,
            verbose=verbose
        )
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])
        if not id_list:
            return [], 0

        runs = []
        expected_runs = len(id_list)  # Proxy
        batch_size = 200
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            handle = safe_entrez_request(
                Entrez.efetch,
                db="sra",
                id=",".join(batch_ids),
                rettype="runinfo",
                retmode="text",
                verbose=verbose
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
                if len(cols) > run_index and cols[run_index].strip():
                    runs.append(cols[run_index].strip())

        # Check XML for SRS and runs using ElementTree
        if verbose:
            with print_lock:
                print(f" - Checking XML for SRS and runs in fallback for {accession}...")
        srs_ids = []
        xml_runs = []
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            handle = safe_entrez_request(
                Entrez.efetch,
                db="sra",
                id=",".join(batch_ids),
                rettype="full",
                retmode="xml",
                verbose=verbose
            )
            xml_text = handle.read()
            handle.close()

            if isinstance(xml_text, bytes):
                xml_text = xml_text.decode("utf-8")

            try:
                root = ET.fromstring(xml_text)
                # Parse XML for SRS and runs
                for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
                    # Extract SRS
                    sample = exp_pkg.find(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
                    if sample is not None and sample.text and sample.text.startswith("SRS"):
                        srs_ids.append(sample.text)
                    # Extract runs from RUN_SET
                    run_set = exp_pkg.find(".//RUN_SET")
                    if run_set is not None:
                        for run in run_set.findall("RUN"):
                            run_acc = run.get("accession")
                            if run_acc and run_acc.startswith(("SRR", "ERR", "DRR")):
                                xml_runs.append(run_acc)
            except ET.ParseError as e:
                with print_lock:
                    print(f"XML parsing failed for {accession}: {e}")
                continue

        if xml_runs:
            if verbose:
                with print_lock:
                    print(f" - Found runs in XML RUN_SET for {accession}.")
            runs.extend(xml_runs)
            expected_runs = max(expected_runs, len(set(xml_runs)))

        if srs_ids:
            if verbose:
                with print_lock:
                    print(f" - Found SRS for {accession}. Fetching runs...")
            for srs in srs_ids:
                handle = safe_entrez_request(
                    Entrez.esearch,
                    db="sra",
                    term=f"{srs}[Sample]",
                    retmax=100000,
                    verbose=verbose
                )
                srs_record = Entrez.read(handle)
                handle.close()
                srs_id_list = srs_record.get("IdList", [])
                if srs_id_list:
                    for j in range(0, len(srs_id_list), batch_size):
                        srs_batch_ids = srs_id_list[j:j + batch_size]
                        handle = safe_entrez_request(
                            Entrez.efetch,
                            db="sra",
                            id=",".join(srs_batch_ids),
                            rettype="runinfo",
                            retmode="text",
                            verbose=verbose
                        )
                        srs_text = handle.read()
                        handle.close()

                        if isinstance(srs_text, bytes):
                            srs_text = srs_text.decode("utf-8")
                        srs_lines = srs_text.strip().split("\n")
                        if len(srs_lines) < 2:
                            continue

                        srs_header = srs_lines[0].split(",")
                        try:
                            srs_run_index = srs_header.index("Run")
                        except ValueError:
                            srs_run_index = 0

                        for line in srs_lines[1:]:
                            cols = line.split(",")
                            if len(cols) > srs_run_index and cols[srs_run_index].strip():
                                runs.append(cols[srs_run_index].strip())

        expected_runs = max(expected_runs, len(set(runs)))
        return list(set(runs)), expected_runs
    except Exception as e:
        with print_lock:
            print(f"Fallback Entrez failed for {accession}: {e}")
        return [], 0

def process_accession(accession, email, api_key, verbose=False):
    """Resolve accession to run accessions (SRR/ERR/DRR) and estimate expected runs."""
    accession = accession.strip()

    # If already a run accession, return it
    if re.match(r"^(SRR|ERR|DRR)\d+$", accession):
        return [accession], 1  # Expected = 1 for direct run accession

    # Primary approach based on prefix
    expected_runs = 0
    if accession.startswith("PRJNA") or accession.startswith("SRP"):
        runs_1, ncbi_expected = get_run_accessions_ncbi(accession, email, api_key, verbose=verbose)
        runs_2 = fetch_runs_ena(accession, verbose=verbose)
        merged_runs = list(set(runs_1 + runs_2))
        expected_runs = ncbi_expected  # Use NCBI's estimate
    elif accession.startswith("GSE"):
        merged_runs, _ = get_run_accessions_geo(accession, email, api_key, verbose=verbose)
        expected_runs = len(merged_runs)  # No reliable estimate, use found runs
    elif (accession.startswith("PRJEB") or
          accession.startswith("ERP") or
          accession.startswith("EGAS")):
        merged_runs = fetch_runs_ena(accession, verbose=verbose)
        expected_runs = len(merged_runs)  # ENA doesn't provide expected count
    else:
        merged_runs = []
        expected_runs = 0

    if merged_runs:
        return merged_runs, expected_runs

    # Fallback with Entrez [ACCN]
    if verbose:
        with print_lock:
            print(f" - No runs found in primary approach for {accession}, trying fallback_entrez_accession...")
    fallback_runs, fallback_expected = fallback_entrez_accession(accession, email, api_key, verbose=verbose)
    expected_runs = max(expected_runs, fallback_expected)
    return fallback_runs, expected_runs

def process_accession_wrapper(args):
    """Wrapper for threading: process one accession and return structured result."""
    accession, email, api_key, verbose = args
    try:
        run_list, expected_runs = process_accession(accession, email, api_key, verbose)
        found_runs = len(run_list)
        status = "Complete" if found_runs >= expected_runs and expected_runs > 0 else "Incomplete" if expected_runs > 0 else "Unknown"
        if verbose:
            with print_lock:
                print(f" - Expected {expected_runs} runs, found {found_runs} runs: {status}")
        return {
            "accession": accession,
            "runs": run_list,
            "expected": expected_runs,
            "found": found_runs,
            "status": status,
            "error": None
        }
    except Exception as e:
        with print_lock:
            print(f"❌ Failed to process {accession}: {e}")
        return {
            "accession": accession,
            "runs": [],
            "expected": 0,
            "found": 0,
            "status": "Failed",
            "error": str(e)
        }

def read_accessions(input_file):
    """Read accessions from input file."""
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

def write_outputs(results, output_file, runs_only_file, fail_log, comparison_log_file):
    """Write all output files."""
    # Sort results by accession to maintain input order
    results.sort(key=lambda x: x["accession"])

    # Collect data
    csv_results = []
    runs_only = []
    failed = []
    comparison = []
    seen_runs = set()

    for result in results:
        accession = result["accession"]
        runs = result["runs"]
        expected = result["expected"]
        found = result["found"]
        status = result["status"]
        error = result["error"]

        comparison.append({
            "accession": accession,
            "expected": expected,
            "found": found,
            "status": status
        })

        if error:
            failed.append((accession, error))
        else:
            for run in runs:
                if run not in seen_runs:
                    csv_results.append((accession, run))
                    runs_only.append(run)
                    seen_runs.add(run)
            if not runs and status != "Failed":
                failed.append((accession, "No run accessions found"))

    # Write CSV
    try:
        with open(output_file, "w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Project Accession", "Run Accession"])
            writer.writerows(csv_results)
        print(f"✅ Results saved to {output_file}")
    except Exception as e:
        print(f"Error writing CSV: {e}")

    # Write runs-only text
    try:
        with open(runs_only_file, "w") as txtfile:
            for run in runs_only:
                txtfile.write(run + "\n")
        print(f"✅ Run accessions only saved to {runs_only_file}")
    except Exception as e:
        print(f"Error writing run accessions only file: {e}")

    # Write failures
    if failed:
        try:
            with open(fail_log, "w") as flog:
                for acc, reason in failed:
                    flog.write(f"{acc}\t{reason}\n")
            print(f"⚠️ Failed accessions logged to {fail_log}")
        except Exception as e:
            print(f"Error writing failure log: {e}")

    # Write comparison log
    try:
        with open(comparison_log_file, "w") as clog:
            clog.write("Accession\tExpected Runs\tFound Runs\tStatus\n")
            for entry in comparison:
                clog.write(f"{entry['accession']}\t{entry['expected']}\t{entry['found']}\t{entry['status']}\n")
        print(f"📊 Run comparison logged to {comparison_log_file}")
    except Exception as e:
        print(f"Error writing comparison log: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Download run accessions for a list of project (or other) accessions, with parallel processing and checkpointing.",
        epilog="Example usage: python accession.py projects.txt --email your_email@example.com --api-key YOUR_API_KEY --threads 4 --verbose"
    )
    parser.add_argument("input_file", help="Input file with one accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="CSV output file name.")
    parser.add_argument("--runs-only", default="run_accessions_only.txt", help="Text file with only run accessions.")
    parser.add_argument("--fail-log", default="failed_accessions.log", help="File to log failed accessions.")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for detailed logging.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing (1 to 10).")
    args = parser.parse_args()

    email = args.email if args.email else input("Enter your email (required for NCBI queries): ").strip()
    api_key = args.api_key if args.api_key else input("Enter your NCBI API key (press enter if none): ").strip()

    if not email:
        print("Error: Email is required for NCBI queries.")
        sys.exit(1)

    # Validate threads
    threads = max(1, min(args.threads, 10))
    if args.threads != threads:
        print(f"Threads adjusted to {threads} (must be between 1 and 10).")

    # Register signal handler for Ctrl+C
    signal.signal(signal.SIGINT, signal_handler)

    # Load checkpoint
    load_checkpoint()

    # Read accessions
    accessions = read_accessions(args.input_file)
    total_accessions = len(accessions)

    # Filter out completed accessions
    completed = set(checkpoint_data["completed"])
    remaining_accessions = [acc for acc in accessions if acc not in completed]
    print(f"Processing {len(remaining_accessions)} of {total_accessions} accessions (skipping {len(completed)} completed).")

    if not remaining_accessions:
        print("All accessions already processed. Writing outputs from checkpoint.")
        write_outputs(checkpoint_data["comparison"], args.output_file, args.runs_only, args.fail_log, "run_comparison.log")
        return

    # Prepare results
    results = checkpoint_data["comparison"].copy()

    # Process accessions
    if threads == 1:
        # Single-threaded with progress bar
        if not args.verbose and tqdm is not None:
            iterator = tqdm(remaining_accessions, desc="Processing accessions", unit="accession")
        else:
            iterator = remaining_accessions

        for accession in iterator:
            if args.verbose:
                with print_lock:
                    print(f"\n🔍 Processing: {accession}")
            result = process_accession_wrapper((accession, email, api_key, args.verbose))
            results.append(result)

            # Update checkpoint
            with checkpoint_lock:
                checkpoint_data["completed"].append(accession)
                checkpoint_data["comparison"].append(result)
                if result["error"]:
                    checkpoint_data["failed"].append((accession, result["error"]))
                else:
                    for run in result["runs"]:
                        if run not in set(checkpoint_data["runs_only"]):
                            checkpoint_data["results"].append((accession, run))
                            checkpoint_data["runs_only"].append(run)
                    if not result["runs"] and result["status"] != "Failed":
                        checkpoint_data["failed"].append((accession, "No run accessions found"))
                save_checkpoint()

            time.sleep(0.5)  # Rate limit
    else:
        # Parallel processing
        tasks = [(acc, email, api_key, args.verbose) for acc in remaining_accessions]
        processed_count = len(completed)

        if not args.verbose:
            print(f"Starting parallel processing with {threads} threads...")

        with ThreadPoolExecutor(max_workers=threads) as executor:
            for result in executor.map(process_accession_wrapper, tasks):
                results.append(result)
                processed_count += 1

                # Update checkpoint
                with checkpoint_lock:
                    checkpoint_data["completed"].append(result["accession"])
                    checkpoint_data["comparison"].append(result)
                    if result["error"]:
                        checkpoint_data["failed"].append((result["accession"], result["error"]))
                    else:
                        for run in result["runs"]:
                            if run not in set(checkpoint_data["runs_only"]):
                                checkpoint_data["results"].append((result["accession"], run))
                                checkpoint_data["runs_only"].append(run)
                        if not result["runs"] and result["status"] != "Failed":
                            checkpoint_data["failed"].append((result["accession"], "No run accessions found"))
                    save_checkpoint()

                # Show progress in non-verbose mode
                if not args.verbose:
                    with print_lock:
                        print(f"Processed {processed_count}/{total_accessions} accessions")

                time.sleep(0.1)  # Slight delay to avoid overwhelming servers

    # Write final outputs
    write_outputs(results, args.output_file, args.runs_only, args.fail_log, "run_comparison.log")

if __name__ == '__main__':
    main()
