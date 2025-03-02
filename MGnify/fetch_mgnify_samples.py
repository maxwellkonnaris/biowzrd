import requests
import pandas as pd
import argparse
import sys
import time
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

LOG_FILE = "fetch_log.txt"

def log_message(message):
    """Append a message to fetch_log.txt (for consistent logging)."""
    with open(LOG_FILE, "a") as log_file:
        log_file.write(message + "\n")

# Define available experiment types
experiment_options = {
    "metagenomic": "Shotgun Metagenomics",
    "16s-rrna-gene-amplicon": "16S rRNA Gene Amplicon (Default)",
    "18s-rrna-gene-amplicon": "18S rRNA Gene Amplicon",
    "its-gene-amplicon": "ITS Gene Amplicon (Fungal)"
}

# A global counter to track how many samples we have processed so far.
processed_count = 0
# A lock to avoid race conditions when threads update processed_count or write logs.
processed_lock = threading.Lock()
# We'll store the time when we start fetching (set in main()).
start_time = None

def fetch_page(url):
    """
    Fetch a single page's worth of sample data from MGnify and parse the JSON
    into our standard metadata dict format.
    """
    global processed_count

    # Perform the GET request
    resp = requests.get(url)
    if resp.status_code != 200:
        log_message(f"❌ Error fetching data from {url}: HTTP {resp.status_code}")
        return []

    data = resp.json()
    samples = data.get("data", [])
    page_results = []
    for sample in samples:
        attributes = sample.get("attributes", {})
        relationships = sample.get("relationships", {})
        study_info = relationships.get("study", {}).get("data", {})
        metadata = {
            "sample_accession": sample.get("id", "N/A"),
            "biome": attributes.get("biome", "N/A"),
            "environment": attributes.get("environment_material", "N/A"),
            "temperature": attributes.get("environment_temperature", "N/A"),
            "salinity": attributes.get("environment_salinity", "N/A"),
            "pH": attributes.get("environment_ph", "N/A"),
            "latitude": attributes.get("latitude", "N/A"),
            "longitude": attributes.get("longitude", "N/A"),
            "collection_date": attributes.get("collection_date", "N/A"),
            "study_accession": study_info.get("id", "N/A"),
            "experiment_type": attributes.get("experiment_type", "N/A"),
        }
        page_results.append(metadata)

    # Update our global progress counter and log if we cross another 1000 boundary
    with processed_lock:
        old_count = processed_count
        processed_count += len(page_results)
        # If, for instance, old_count=2999 and new_count=3003, we've crossed 3000
        # So let's check if we've crossed a multiple of 1000 in this batch
        check_range = range(old_count // 1000 + 1, processed_count // 1000 + 1)
        if len(check_range) > 0:  # We crossed at least one multiple of 1000
            elapsed_time = time.time() - start_time
            log_message(f"Processed {processed_count} samples... (Elapsed Time: {elapsed_time:.2f} seconds)")

    return page_results

def main():
    global start_time
    parser = argparse.ArgumentParser(
        description="Fetch bulk metadata from MGnify and save it as a CSV file."
    )
    parser.add_argument(
        "--experiment",
        type=str,
        choices=experiment_options.keys(),
        default="16s-rrna-gene-amplicon",
        help=("Specify the experiment type. Available options: "
              "metagenomic, 16s-rrna-gene-amplicon (default), "
              "18s-rrna-gene-amplicon, its-gene-amplicon.")
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of threads to use in parallel fetching."
    )
    args = parser.parse_args()
    experiment_type = args.experiment
    max_workers = args.threads

    # 1) Clear previous logs and write a "start" message
    with open(LOG_FILE, "w") as log_file:
        log_file.write("=== Bulk Fetch Script Execution Started ===\n")
    log_message(f"Experiment type: {experiment_options[experiment_type]}")
    log_message(f"Threads: {max_workers}")

    print("\n🔹 **MGnify Bulk Metadata Fetcher** 🔹")
    print(f"✅ Fetching bulk metadata for: {experiment_options[experiment_type]}")
    print(f"📌 Progress is logged to `{LOG_FILE}`.")

    # 2) We do an initial request just to discover the total count (and maybe grab page=1).
    base_url = "https://www.ebi.ac.uk/metagenomics/api/latest/samples"
    page_size = 1000  # or 2000, 5000, etc. Fine-tune as needed.
    # Page 1 request
    first_page_url = f"{base_url}?experiment-type={experiment_type}&page_size={page_size}&page=1"
    resp = requests.get(first_page_url)
    if resp.status_code != 200:
        error_message = f"❌ Error fetching first page: {resp.status_code}"
        print(error_message)
        log_message(error_message)
        sys.exit(1)

    data = resp.json()
    meta_pagination = data.get("meta", {}).get("pagination", {})
    total_count = meta_pagination.get("count")
    if not total_count:
        print("Could not determine 'count' from API response; aborting.")
        log_message("Could not determine 'count' from API. Aborting.")
        sys.exit(1)

    num_pages = math.ceil(total_count / page_size)
    log_message(f"Total samples reported: {total_count}")
    log_message(f"Total pages needed: {num_pages}")

    # 3) Start concurrency to fetch all pages, including page=1
    all_results = []
    all_page_urls = []
    for page_idx in range(1, num_pages + 1):
        url = f"{base_url}?experiment-type={experiment_type}&page_size={page_size}&page={page_idx}"
        all_page_urls.append(url)

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
    print(f"\n✅ Retrieved metadata for {len(all_results)} samples in {total_time:.2f} seconds.")
    log_message(f"✅ Completed fetching metadata for {len(all_results)} samples in {total_time:.2f} seconds.")

    # 4) Convert metadata to a DataFrame and save to CSV
    df = pd.DataFrame(all_results)
    csv_filename = f"mgnify_bulk_metadata_{experiment_type.replace('-', '_')}.csv"
    df.to_csv(csv_filename, index=False)

    print(f"\n✅ **Metadata saved to '{csv_filename}'.** 🎉")
    log_message(f"✅ Metadata saved to '{csv_filename}'.")
    print("\n🎯 Bulk fetching completed. Check `fetch_log.txt` for details.")
    log_message("🎯 Bulk fetching completed successfully.")


if __name__ == "__main__":
    main()
