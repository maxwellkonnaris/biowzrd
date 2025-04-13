#!/usr/bin/env python3
import argparse
import csv
import sys
import time
import requests
import re
import xml.etree.ElementTree as ET

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
    print("Warning: tqdm not installed. Progress bar will be disabled. Install with 'pip install tqdm' for progress display.")

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
                    print(f"Bad request (400) for {url}: {e}. Retrying in {wait_time}s...")
            else:
                if verbose:
                    print(f"Request error: {e}. Retrying in {wait_time}s...")
            time.sleep(wait_time)
        except requests.exceptions.RequestException as e:
            wait_time = 2 ** attempt
            if verbose:
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
                print(f"ENA fetch failed for {accession}: Invalid accession or no runs available.")
        else:
            if verbose:
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
                print(f"XML parsing failed for {accession}: {e}")
                continue

        if xml_runs:
            if verbose:
                print(f" - Found runs in XML RUN_SET for {accession}: {', '.join(xml_runs)}")
            runs.extend(xml_runs)
            # Update expected runs based on XML RUN_SET
            expected_runs = max(expected_runs, len(set(xml_runs)))

        # Step 4: Query SRA for runs linked to SRS
        if srs_ids:
            if verbose:
                print(f" - Found SRS for {accession}: {', '.join(srs_ids)}. Fetching runs...")
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
        print(f"Entrez failed for GEO accession {accession}: {e}")

    if not runs and SRAweb is not None:
        db = SRAweb()
        for attempt in range(2):
            try:
                if verbose:
                    print(f"Falling back to pysradb for GEO accession {accession} (attempt {attempt + 1})")
                df = db.sra_metadata(geo=accession, detailed=True)
                if not df.empty and "run_accession" in df.columns:
                    runs.extend(df["run_accession"].dropna().unique().tolist())
                break
            except Exception as e:
                wait_time = 2 ** attempt
                if verbose:
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
                print(f"XML parsing failed for {accession}: {e}")
                continue

        if xml_runs:
            if verbose:
                print(f" - Found runs in XML RUN_SET for {accession}: {', '.join(xml_runs)}")
            runs.extend(xml_runs)
            expected_runs = max(expected_runs, len(set(xml_runs)))

        if srs_ids:
            if verbose:
                print(f" - Found SRS for {accession}: {', '.join(srs_ids)}. Fetching runs...")
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
        print(f" - No runs found in primary approach for {accession}, trying fallback_entrez_accession...")
    fallback_runs, fallback_expected = fallback_entrez_accession(accession, email, api_key, verbose=verbose)
    expected_runs = max(expected_runs, fallback_expected)
    return fallback_runs, expected_runs

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

def main():
    parser = argparse.ArgumentParser(
        description="Download run accessions for a list of project (or other) accessions sequentially.",
        epilog="Example usage: python download_runs.py projects.txt --email your_email@example.com --api-key YOUR_API_KEY --verbose"
    )
    parser.add_argument("input_file", help="Input file with one accession per line.")
    parser.add_argument("-o", "--output_file", default="run_accessions.csv", help="CSV output file name.")
    parser.add_argument("--runs-only", default="run_accessions_only.txt", help="Text file with only run accessions.")
    parser.add_argument("--fail-log", default="failed_accessions.log", help="File to log failed accessions.")
    parser.add_argument("--email", help="Your email address (required for NCBI queries).")
    parser.add_argument("--api-key", help="Your NCBI API key (optional).")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for detailed logging.")
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
    comparison_log = []

    # Process accessions with progress bar in non-verbose mode
    if not args.verbose and tqdm is not None:
        iterator = tqdm(accessions, desc="Processing accessions", unit="accession")
    else:
        iterator = accessions

    for accession in iterator:
        if args.verbose:
            print(f"\n🔍 Processing: {accession}")
        try:
            run_list, expected_runs = process_accession(accession, email, api_key, verbose=args.verbose)
            found_runs = len(run_list)
            status = "Complete" if found_runs >= expected_runs and expected_runs > 0 else "Incomplete" if expected_runs > 0 else "Unknown"
            comparison_log.append({
                "accession": accession,
                "expected": expected_runs,
                "found": found_runs,
                "runs": run_list,
                "status": status
            })
            if args.verbose:
                print(f" - Expected {expected_runs} runs, found {found_runs} runs: {status}")
                if run_list:
                    print(f" - Runs: {', '.join(run_list)}")

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
            comparison_log.append({
                "accession": accession,
                "expected": 0,
                "found": 0,
                "runs": [],
                "status": "Failed"
            })
        time.sleep(0.5)

    # Write CSV
    try:
        with open(args.output_file, "w", newline="", encoding="utf-8-sig") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Project Accession", "Run Accession"])
            writer.writerows(results)
        print(f"\n✅ Results saved to {args.output_file}")
    except Exception as e:
        print(f"Error writing CSV: {e}")

    # Write runs-only text
    try:
        with open(args.runs_only, "w") as txtfile:
            for run in runs_only:
                txtfile.write(run + "\n")
        print(f"✅ Run accessions only saved to {args.runs_only}")
    except Exception as e:
        print(f"Error writing run accessions only file: {e}")

    # Write failures
    if failed:
        try:
            with open(args.fail_log, "w") as flog:
                for acc, reason in failed:
                    flog.write(f"{acc}\t{reason}\n")
            print(f"⚠️ Failed accessions logged to {args.fail_log}")
        except Exception as e:
            print(f"Error writing failure log: {e}")

    # Write comparison log
    try:
        with open("run_comparison.log", "w") as clog:
            clog.write("Accession\tExpected Runs\tFound Runs\tStatus\tRun Accessions\n")
            for entry in comparison_log:
                runs_str = ",".join(entry["runs"]) if entry["runs"] else "None"
                clog.write(f"{entry['accession']}\t{entry['expected']}\t{entry['found']}\t{entry['status']}\t{runs_str}\n")
        print(f"📊 Run comparison logged to run_comparison.log")
    except Exception as e:
        print(f"Error writing comparison log: {e}")

if __name__ == '__main__':
    main()
