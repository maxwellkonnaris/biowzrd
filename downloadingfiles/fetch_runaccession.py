#!/usr/bin/env python3
import argparse
import csv
import sys
import time
import requests
import re
import xml.etree.ElementTree as ET
import json
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
import shutil

# Store current working directory
CWD = os.getcwd()
print(f"DEBUG: Current working directory is {CWD}", flush=True)

# Check if CWD is writable
try:
    with open(os.path.join(CWD, "test_write.txt"), "w") as f:
        f.write("test")
    os.remove(os.path.join(CWD, "test_write.txt"))
    print("DEBUG: CWD is writable", flush=True)
except Exception as e:
    print(f"Error: CWD {CWD} is not writable: {e}", flush=True)
    sys.exit(1)

# Force output flushing for SLURM
sys.stdout.flush()
print("DEBUG: Script initializing", flush=True)

# Check if Biopython is installed.
try:
    from Bio import Entrez
except ImportError:
    print("Error: Biopython is required but not installed. Please install it (e.g., pip install biopython) and try again.", flush=True)
    sys.exit(1)

# Attempt to import pysradb for fallback on GEO accessions.
try:
    from pysradb import SRAweb
except ImportError:
    SRAweb = None
    print("DEBUG: pysradb not installed, GEO fallback disabled", flush=True)

# Attempt to import tqdm for progress bar.
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None
    print("Warning: tqdm not installed. Progress bar disabled. Install with 'pip install tqdm' for progress display.", flush=True)

def safe_entrez_request(func, *args, max_retries=3, verbose=False, **kwargs):
    """Helper to gracefully retry Entrez calls on HTTP 429 or transient failures."""
    for attempt in range(max_retries):
        try:
            handle = func(*args, **kwargs)
            return handle
        except Exception as e:
            if "HTTP Error 429" in str(e):
                wait_time = 5 * 2 ** attempt
                if verbose:
                    print(f"Entrez 429: Retrying in {wait_time}s...", flush=True)
                time.sleep(wait_time)
            else:
                raise e
    raise RuntimeError(f"Entrez request failed after {max_retries} retries.")

def safe_requests_get(url, params=None, max_retries=3, timeout=30, verbose=False):
    """Helper to retry requests.get on transient failures."""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.exceptions.HTTPError as e:
            wait_time = 5 * 2 ** attempt
            if e.response.status_code == 400:
                if verbose:
                    print(f"Bad request (400) for {url}: {e}. Retrying in {wait_time}s...", flush=True)
            else:
                if verbose:
                    print(f"Request error: {e}. Retrying in {wait_time}s...", flush=True)
            time.sleep(wait_time)
        except requests.exceptions.RequestException as e:
            wait_time = 5 * 2 ** attempt
            if verbose:
                print(f"Request error: {e}. Retrying in {wait_time}s...", flush=True)
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
                print(f"ENA fetch failed for {accession}: Invalid accession or no runs available.", flush=True)
        else:
            if verbose:
                print(f"ENA fetch failed for {accession}: {e}", flush=True)
        return []

def get_run_accessions_ncbi(accession, email, api_key, verbose=False):
    """Fetch run accessions from NCBI SRA for BioProject (PRJNA/SRP), including SRS and XML RUN_SET."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        if verbose:
            print(f" - Searching SRA for BioProject {accession}...", flush=True)
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
                print(f" - No SRA experiments found for {accession}.", flush=True)
            return [], 0

        expected_runs = len(id_list)
        runs = []
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

        if verbose:
            print(f" - Checking XML for SRS and runs in {accession}...", flush=True)
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
                for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
                    sample = exp_pkg.find(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
                    if sample is not None and sample.text and sample.text.startswith("SRS"):
                        srs_ids.append(sample.text)
                    run_set = exp_pkg.find(".//RUN_SET")
                    if run_set is not None:
                        for run in run_set.findall("RUN"):
                            run_acc = run.get("accession")
                            if run_acc and run_acc.startswith(("SRR", "ERR", "DRR")):
                                xml_runs.append(run_acc)
            except ET.ParseError as e:
                print(f"XML parsing failed for {accession}: {e}", flush=True)
                continue

        if xml_runs:
            if verbose:
                print(f" - Found runs in XML RUN_SET for {accession}.", flush=True)
            runs.extend(xml_runs)
            expected_runs = max(expected_runs, len(set(xml_runs)))

        if srs_ids:
            if verbose:
                print(f" - Found SRS for {accession}. Fetching runs...", flush=True)
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
        print(f"Error processing NCBI accession {accession}: {e}", flush=True)
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
        print(f"Entrez failed for GEO accession {accession}: {e}", flush=True)

    if not runs and SRAweb is not None:
        db = SRAweb()
        for attempt in range(2):
            try:
                if verbose:
                    print(f"Falling back to pysradb for GEO accession {accession} (attempt {attempt + 1})", flush=True)
                df = db.sra_metadata(geo=accession, detailed=True)
                if not df.empty and "run_accession" in df.columns:
                    runs.extend(df["run_accession"].dropna().unique().tolist())
                break
            except Exception as e:
                wait_time = 5 * 2 ** attempt
                if verbose:
                    print(f"pysradb error: {e}. Retrying in {wait_time}s...", flush=True)
                time.sleep(wait_time)

    return list(set(runs)), 0

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
        expected_runs = len(id_list)
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

        if verbose:
            print(f" - Checking XML for SRS and runs in fallback for {accession}...", flush=True)
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
                for exp_pkg in root.findall(".//EXPERIMENT_PACKAGE"):
                    sample = exp_pkg.find(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
                    if sample is not None and sample.text and sample.text.startswith("SRS"):
                        srs_ids.append(sample.text)
                    run_set = exp_pkg.find(".//RUN_SET")
                    if run_set is not None:
                        for run in run_set.findall("RUN"):
                            run_acc = run.get("accession")
                            if run_acc and run_acc.startswith(("SRR", "ERR", "DRR")):
                                xml_runs.append(run_acc)
            except ET.ParseError as e:
                print(f"XML parsing failed for {accession}: {e}", flush=True)
                continue

        if xml_runs:
            if verbose:
                print(f" - Found runs in XML RUN_SET for {accession}.", flush=True)
            runs.extend(xml_runs)
            expected_runs = max(expected_runs, len(set(xml_runs)))

        if srs_ids:
            if verbose:
                print(f" - Found SRS for {accession}. Fetching runs...", flush=True)
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
        print(f"Fallback Entrez failed for {accession}: {e}", flush=True)
        return [], 0

def process_accession(accession, email, api_key, verbose=False):
    """Resolve accession to run accessions (SRR/ERR/DRR) and estimate expected runs."""
    accession = accession.strip()
    if verbose:
        print(f"DEBUG: Processing accession {accession}", flush=True)

    if re.match(r"^(SRR|ERR|DRR)\d+$", accession):
        return [accession], 1

    expected_runs = 0
    if accession.startswith("PRJNA") or accession.startswith("SRP"):
        runs_1, ncbi_expected = get_run_accessions_ncbi(accession, email, api_key, verbose=verbose)
        runs_2 = fetch_runs_ena(accession, verbose=verbose)
        merged_runs = list(set(runs_1 + runs_2))
        expected_runs = ncbi_expected
    elif accession.startswith("GSE"):
        merged_runs, _ = get_run_accessions_geo(accession, email, api_key, verbose=verbose)
        expected_runs = len(merged_runs)
    elif (accession.startswith("PRJEB") or
          accession.startswith("ERP") or
          accession.startswith("EGAS")):
        merged_runs = fetch_runs_ena(accession, verbose=verbose)
        expected_runs = len(merged_runs)
    else:
        merged_runs = []
        expected_runs = 0

    if merged_runs:
        return merged_runs, expected_runs

    if verbose:
        print(f" - No runs found in primary approach for {accession}, trying fallback_entrez_accession...", flush=True)
    fallback_runs, fallback_expected = fallback_entrez_accession(accession, email, api_key, verbose=verbose)
    expected_runs = max(expected_runs, fallback_expected)
    return fallback_runs, expected_runs

def process_accession_wrapper(args):
    """Wrapper for threading: process one accession and save to temp file."""
    accession, email, api_key, verbose, temp_dir = args
    try:
        run_list, expected_runs = process_accession(accession, email, api_key, verbose)
        found_runs = len(run_list)
        status = (
            "Complete"
            if found_runs >= expected_runs and expected_runs > 0
            else "Incomplete" if expected_runs > 0
            else "Unknown"
        )
        result = {
            "accession": accession,
            "runs": run_list,
            "expected": expected_runs,
            "found": found_runs,
            "status": status,
            "error": None
        }
        # Save to temp file
        temp_file = os.path.join(temp_dir, f"tmp_{accession}.json")
        try:
            with open(temp_file, "w") as f:
                json.dump(result, f, indent=2)
            if verbose:
                print(f"DEBUG: Saved temp result to {temp_file}", flush=True)
        except Exception as e:
            print(f"Error saving temp file {temp_file}: {e}", flush=True)
        if verbose:
            print(f" - Expected {expected_runs} runs, found {found_runs} runs: {status}", flush=True)
        return result
    except Exception as e:
        print(f"❌ Failed to process {accession}: {e}", flush=True)
        result = {
            "accession": accession,
            "runs": [],
            "expected": 0,
            "found": 0,
            "status": "Failed",
            "error": str(e)
        }
        # Save to temp file even on failure
        temp_file = os.path.join(temp_dir, f"tmp_{accession}.json")
        try:
            with open(temp_file, "w") as f:
                json.dump(result, f, indent=2)
            if verbose:
                print(f"DEBUG: Saved temp result (error) to {temp_file}", flush=True)
        except Exception as e:
            print(f"Error saving temp file {temp_file}: {e}", flush=True)
        return result

def read_accessions(input_file):
    """Read accessions from input file."""
    input_path = os.path.join(CWD, input_file)
    print(f"DEBUG: Reading input file {input_path}", flush=True)
    accessions = []
    try:
        with open(input_path, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    accessions.append(line)
        print(f"DEBUG: Read {len(accessions)} accessions from {input_path}", flush=True)
    except Exception as e:
        print(f"Error reading {input_path}: {e}", flush=True)
        sys.exit(1)
    return accessions

def merge_results(results, temp_dir):
    """Merge results from memory and temp files, ensuring no data loss."""
    print(f"DEBUG: Merging results from {temp_dir}", flush=True)
    merged_results = results.copy()
    try:
        for temp_file in os.listdir(temp_dir):
            if temp_file.startswith("tmp_") and temp_file.endswith(".json"):
                temp_path = os.path.join(temp_dir, temp_file)
                try:
                    with open(temp_path, "r") as f:
                        temp_result = json.load(f)
                        if not any(r["accession"] == temp_result["accession"] for r in merged_results):
                            merged_results.append(temp_result)
                            print(f"DEBUG: Added temp result from {temp_path}", flush=True)
                except Exception as e:
                    print(f"Error reading temp file {temp_path}: {e}", flush=True)
    except Exception as e:
        print(f"Error accessing temp directory {temp_dir}: {e}", flush=True)
    print(f"DEBUG: Merged {len(merged_results)} results", flush=True)
    return merged_results

def write_outputs(results, output_file, runs_only_file, fail_log, comparison_log_file):
    """Write all output files."""
    output_path = os.path.join(CWD, output_file)
    runs_only_path = os.path.join(CWD, runs_only_file)
    fail_log_path = os.path.join(CWD, fail_log)
    comparison_log_path = os.path.join(CWD, comparison_log_file)
    
    print(f"DEBUG: Starting write_outputs to {CWD} with {len(results)} results", flush=True)
    results.sort(key=lambda x: x["accession"])

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

    try:
        print(f"DEBUG: Writing CSV to {output_path}", flush=True)
        with open(output_path, "w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Project Accession", "Run Accession"])
            writer.writerows(csv_results)
        print(f"✅ Results saved to {output_path}", flush=True)
    except Exception as e:
        print(f"Error writing CSV to {output_path}: {e}", flush=True)

    try:
        print(f"DEBUG: Writing runs-only to {runs_only_path}", flush=True)
        with open(runs_only_path, "w") as txtfile:
            for run in runs_only:
                txtfile.write(run + "\n")
        print(f"✅ Run accessions only saved to {runs_only_path}", flush=True)
    except Exception as e:
        print(f"Error writing run accessions only file to {runs_only_path}: {e}", flush=True)

    if failed:
        try:
            print(f"DEBUG: Writing failure log to {fail_log_path}", flush=True)
            with open(fail_log_path, "w") as flog:
                for acc, reason in failed:
                    flog.write(f"{acc}\t{reason}\n")
            print(f"⚠️ Failed accessions logged to {fail_log_path}", flush=True)
        except Exception as e:
            print(f"Error writing failure log to {fail_log_path}: {e}", flush=True)

    try:
        print(f"DEBUG: Writing comparison log to {comparison_log_path}", flush=True)
        with open(comparison_log_path, "w") as clog:
            clog.write("Accession\tExpected Runs\tFound Runs\tStatus\n")
            for entry in comparison:
                clog.write(f"{entry['accession']}\t{entry['expected']}\t{entry['found']}\t{entry['status']}\n")
        print(f"📊 Run comparison logged to {comparison_log_path}", flush=True)
    except Exception as e:
        print(f"Error writing comparison log to {comparison_log_path}: {e}", flush=True)

def main():
    print("DEBUG: Entering main()", flush=True)
    parser = argparse.ArgumentParser(
        description="Download run accessions for a list of project (or other) accessions, "
                    "with parallel processing.",
        epilog="Example usage: python accession.py studies.txt --email your_email@example.com "
               "--api-key YOUR_API_KEY --threads 8 --verbose"
    )
    parser.add_argument("input_file", help="Input file with one accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="CSV output file name.")
    parser.add_argument("--runs-only", default="run_accessions_only.txt", help="Text file with only run accessions.")
    parser.add_argument("--fail-log", default="failed_accessions.log", help="File to log failed accessions.")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for detailed logging.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for parallel processing (1 to 10).")
    parser.add_argument("--slurm", action="store_true",
                        help="Submit as a SLURM batch job instead of running interactively.")
    args = parser.parse_args()
    print(f"DEBUG: Parsed arguments: {vars(args)}", flush=True)

    # Handle SLURM batch submission
    if args.slurm:
        import subprocess
        print("DEBUG: Preparing SLURM submission", flush=True)
        slurm_script = f"""#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --output={os.path.join(CWD, 'fetch-%j.out')}
#SBATCH --error={os.path.join(CWD, 'fetch-%j.err')}
#SBATCH --chdir={CWD}

echo "DEBUG: SLURM job starting"
echo "DEBUG: Current working directory: $PWD"
echo "DEBUG: Python version: $(python --version)"
echo "DEBUG: Python path: $(which python)"
python {os.path.abspath(__file__)} {args.input_file} -o {args.output_file} --runs-only {args.runs_only} --fail-log {args.fail_log} --email {args.email or ''} --api-key {args.api_key or ''} {'--verbose' if args.verbose else ''} --threads {args.threads}
"""
        slurm_script_path = os.path.join(CWD, "accession_slurm.sh")
        try:
            with open(slurm_script_path, "w") as f:
                f.write(slurm_script)
            print(f"DEBUG: Wrote SLURM script to {slurm_script_path}", flush=True)
            result = subprocess.run(["sbatch", slurm_script_path], capture_output=True, text=True, check=True)
            print(f"Submitted SLURM job: {result.stdout.strip()}", flush=True)
            sys.exit(0)
        except subprocess.CalledProcessError as e:
            print(f"Error submitting SLURM job: {e.stderr}", flush=True)
            sys.exit(1)

    email = args.email if args.email else input("Enter your email (required for NCBI queries): ").strip()
    api_key = args.api_key if args.api_key else input("Enter your NCBI API key (press enter if none): ").strip()

    if not email:
        print("Error: Email is required for NCBI queries.", flush=True)
        sys.exit(1)

    # Validate threads
    threads = max(1, min(args.threads, 10))
    if args.threads != threads:
        print(f"Threads adjusted to {threads} (must be between 1 and 10).", flush=True)

    # Read accessions
    accessions = read_accessions(args.input_file)
    total_accessions = len(accessions)
    print(f"Processing {total_accessions} accessions.", flush=True)

    # Create temp directory in CWD
    temp_dir = os.path.join(CWD, "temp_accession")
    try:
        os.makedirs(temp_dir, exist_ok=True)
        print(f"DEBUG: Created temp directory {temp_dir}", flush=True)
    except Exception as e:
        print(f"Error creating temp directory {temp_dir}: {e}", flush=True)
        sys.exit(1)

    # Prepare results
    results = []

    # Process accessions
    print(f"DEBUG: Starting processing with {threads} threads", flush=True)
    tasks = [(acc, email, api_key, args.verbose, temp_dir) for acc in accessions]
    if threads == 1:
        # Single-threaded
        print("DEBUG: Running in single-threaded mode", flush=True)
        iterator = tqdm(accessions, desc="Processing accessions", unit="accession") if tqdm and not args.verbose else accessions
        for accession in iterator:
            if args.verbose:
                print(f"\n🔍 Processing: {accession}", flush=True)
            result = process_accession_wrapper((accession, email, api_key, args.verbose, temp_dir))
            results.append(result)
            time.sleep(1.0)
    else:
        # Parallel processing
        print(f"DEBUG: Running in parallel mode with {threads} threads", flush=True)
        pbar = tqdm(total=len(accessions), desc="Processing accessions", unit="accession") if tqdm and not args.verbose else None
        try:
            with ThreadPoolExecutor(max_workers=threads) as executor:
                for i, result in enumerate(executor.map(process_accession_wrapper, tasks)):
                    results.append(result)
                    print(f"DEBUG: Thread completed for accession {result['accession']}", flush=True)
                    if pbar:
                        pbar.update(1)
                    time.sleep(1.0)  # Increased delay
        finally:
            if pbar:
                pbar.close()
        print("DEBUG: All threads completed", flush=True)

    # Merge results from memory and temp files
    results = merge_results(results, temp_dir)

    # Write final outputs
    print("DEBUG: Writing final outputs", flush=True)
    write_outputs(results, args.output_file, args.runs_only, args.fail_log, "run_comparison.log")

    # Clean up temp directory
    print(f"DEBUG: Cleaning up temp directory {temp_dir}", flush=True)
    try:
        shutil.rmtree(temp_dir, ignore_errors=True)
        print(f"DEBUG: Temp directory {temp_dir} removed", flush=True)
    except Exception as e:
        print(f"Error cleaning up temp directory {temp_dir}: {e}", flush=True)
    print("DEBUG: Script completed", flush=True)

if __name__ == "__main__":
    print("DEBUG: Script starting", flush=True)
    main()
