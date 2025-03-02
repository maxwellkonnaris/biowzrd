import requests
import pandas as pd
import argparse
import sys
import time
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

LOG_FILE = "fetch_log.txt"

# A dictionary mapping user-friendly experiment keys to the actual param used in `?experiment-type=`
experiment_options = {
    "metagenomic": "metagenomic",
    "16s-rrna-gene-amplicon": "16s-rrna-gene-amplicon",
    "18s-rrna-gene-amplicon": "18s-rrna-gene-amplicon",
    "its-gene-amplicon": "its-gene-amplicon"
}

# Global counters / locks
processed_count = 0
processed_lock = threading.Lock()
start_time = None  # Set at runtime in main()

def log_message(message):
    """Append a message to fetch_log.txt (for consistent logging)."""
    with open(LOG_FILE, "a") as log_file:
        log_file.write(message + "\n")

def parse_sample(sample_json):
    """
    Extract ALL possible metadata from the 'attributes' dict plus 
    everything in the 'sample-metadata' array, plus a couple items from 'relationships'.
    Return a single dict for the DataFrame.
    """
    row = {}

    # Basic ID
    row["sample_id"] = sample_json.get("id", "N/A")

    # Copy everything from "attributes" except "sample-metadata"
    attributes = sample_json.get("attributes", {}).copy()
    sample_md = attributes.pop("sample-metadata", [])

    for attr_key, attr_val in attributes.items():
        row[attr_key] = attr_val

    # Expand "sample-metadata" into columns:
    for md_item in sample_md:
        k = md_item.get("key")
        v = md_item.get("value")
        if k:
            row[k] = v

    # Grab some relationship info (like studies, biome):
    relationships = sample_json.get("relationships", {})
    studies_list = relationships.get("studies", {}).get("data", [])
    if studies_list:
        row["study_ids"] = ";".join([s.get("id", "") for s in studies_list])
    else:
        row["study_ids"] = None

    biome_data = relationships.get("biome", {}).get("data", {})
    row["biome_id"] = biome_data.get("id", None)

    return row

def fetch_page(url):
    """
    Fetch and parse a single page from the MGnify /samples endpoint.
    Returns a list of dicts (one per sample).
    """
    global processed_count

    resp = requests.get(url)
    if resp.status_code != 200:
        log_message(f"❌ Error fetching data from {url}: HTTP {resp.status_code}")
        return []

    data = resp.json()
    samples = data.get("data", [])
    page_results = []

    for sample_json in samples:
        row_dict = parse_sample(sample_json)
        page_results.append(row_dict)

    # Update global progress
    with processed_lock:
        old_count = processed_count
        processed_count += len(page_results)

        # For example, log progress every 1,000 samples
        crossed = range(old_count // 1000 + 1, processed_count // 1000 + 1)
        if len(list(crossed)) > 0:
            elapsed_time = time.time() - start_time
            log_message(f"Processed {processed_count} samples... "
                        f"(Elapsed: {elapsed_time:.2f} seconds)")

    return page_results

def main():
    global start_time

    parser = argparse.ArgumentParser(
        description="Fetch sample metadata from MGnify for a given experiment type, using concurrency."
    )
    parser.add_argument(
        "--experiment",
        type=str,
        choices=experiment_options.keys(),
        default="16s-rrna-gene-amplicon",
        help=("Experiment type. One of: metagenomic, 16s-rrna-gene-amplicon, "
              "18s-rrna-gene-amplicon, its-gene-amplicon.")
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use for parallel fetching."
    )
    parser.add_argument(
        "--page_size",
        type=int,
        default=1000,
        help="Number of samples per page (in the API request)."
    )

    args = parser.parse_args()
    experiment_choice = args.experiment
    max_workers = args.threads
    page_size = args.page_size

    # Convert the CLI experiment arg to the actual query param
    experiment_param = experiment_options[experiment_choice]

    # 1) Wipe log and record start
    with open(LOG_FILE, "w") as log_file:
        log_file.write("=== MGnify Fetch Script Execution Started ===\n")
    log_message(f"Experiment type: {experiment_choice} ({experiment_param})")
    log_message(f"Threads: {max_workers}")
    log_message(f"Page size: {page_size}")

    print("\n🔹 **MGnify Bulk Metadata Fetcher** 🔹")
    print(f"Experiment: {experiment_choice} ({experiment_param})")
    print(f"Using up to {max_workers} concurrent threads.")
    print(f"Page size: {page_size}")
    print(f"Progress will be logged to `{LOG_FILE}`.\n")

    # 2) Start with page=1 and figure out how many total pages are needed
    base_url = "https://www.ebi.ac.uk/metagenomics/api/v1/samples"
    first_url = f"{base_url}?experiment-type={experiment_param}&page=1&page_size={page_size}"
    first_resp = requests.get(first_url)
    if first_resp.status_code != 200:
        error_message = (f"❌ Error fetching first page: {first_url} "
                         f"(HTTP {first_resp.status_code})")
        print(error_message)
        log_message(error_message)
        sys.exit(1)

    first_data = first_resp.json()
    links = first_data.get("links", {})
    last_link = links.get("last")
    if not last_link:
        error_message = ("Could not find 'links[\"last\"]' in the response. Aborting. "
                         "Ensure 'experiment-type' is correct or the endpoint isn't empty.")
        print(error_message)
        log_message(error_message)
        sys.exit(1)

    # last_link might look like:
    #   https://www.ebi.ac.uk/metagenomics/api/v1/samples?experiment-type=16s-rrna-gene-amplicon&page=1234&page_size=1000
    try:
        # Easiest parse: split on "&page=" or parse the query properly
        # We'll do a simple split on "page=" for demonstration:
        after_page_equals = last_link.split("page=")[1]
        last_page_num = int(after_page_equals.split("&")[0])
    except Exception as e:
        error_message = (f"Could not parse the page number from last link: {last_link}\n{e}")
        print(error_message)
        log_message(error_message)
        sys.exit(1)

    # 3) Build the list of all page URLs
    all_page_urls = []
    for page_idx in range(1, last_page_num + 1):
        url = (f"{base_url}?experiment-type={experiment_param}"
               f"&page={page_idx}&page_size={page_size}")
        all_page_urls.append(url)

    total_pages = len(all_page_urls)
    log_message(f"Discovered {total_pages} pages (1..{last_page_num}) for {experiment_choice}.")

    # 4) Fetch them concurrently
    all_results = []
    start_time = time.time()
    log_message("Starting parallel fetch of all pages...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_url = {executor.submit(fetch_page, url): url for url in all_page_urls}
        for future in as_completed(future_to_url):
            url = future_to_url[future]
            try:
                page_data = future.result()
                all_results.extend(page_data)
            except Exception as exc:
                log_message(f"❌ Exception fetching page {url}: {exc}")

    total_time = time.time() - start_time
    final_count = len(all_results)

    print(f"\n✅ Retrieved metadata for {final_count} samples in {total_time:.2f} seconds.")
    log_message(f"✅ Completed fetching metadata for {final_count} samples in {total_time:.2f} seconds.")

    # 5) Convert to DataFrame and save
    df = pd.DataFrame(all_results)
    csv_filename = f"mgnify_samples_{experiment_choice}.csv"
    df.to_csv(csv_filename, index=False)

    print(f"✅ **All sample metadata saved to '{csv_filename}'.** 🎉\n")
    log_message(f"✅ Metadata saved to '{csv_filename}'.")
    print("🎯 Bulk fetching completed. Check `fetch_log.txt` for details.")
    log_message("🎯 Bulk fetching completed successfully.")

if __name__ == "__main__":
    main()
