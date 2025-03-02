import requests
import pandas as pd
import argparse
import sys
import time
import concurrent.futures

LOG_FILE = "fetch_log.txt"

def log_message(message):
    """Write log messages to a file."""
    with open(LOG_FILE, "a") as log_file:
        log_file.write(message + "\n")

# Define available experiment types
experiment_options = {
    "metagenomic": "Shotgun Metagenomics",
    "16s-rrna-gene-amplicon": "16S rRNA Gene Amplicon (Default)",
    "18s-rrna-gene-amplicon": "18S rRNA Gene Amplicon",
    "its-gene-amplicon": "ITS Gene Amplicon (Fungal)"
}

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Fetch metadata from MGnify based on experiment type and save it as a CSV file."
)

parser.add_argument(
    "--experiment",
    type=str,
    choices=experiment_options.keys(),
    default="16s-rrna-gene-amplicon",
    help="Specify the experiment type. Available options: metagenomic, 16s-rrna-gene-amplicon (default), 18s-rrna-gene-amplicon, its-gene-amplicon."
)

args = parser.parse_args()
experiment_type = args.experiment  # Get the selected experiment type

# Clear previous logs
with open(LOG_FILE, "w") as log_file:
    log_file.write("=== Fetch Script Execution Started ===\n")

print("\n🔹 **MGnify Metadata Fetcher** 🔹")
print(f"✅ Fetching metadata for: {experiment_options[experiment_type]}\n")
print(f"📌 Fetching progress is logged to `{LOG_FILE}`.")

def fetch_samples(experiment_type):
    """Fetch ALL available samples from MGnify and log progress."""
    base_url = "https://www.ebi.ac.uk/metagenomics/api/latest/samples"
    sample_accessions = []
    next_page = f"{base_url}?experiment-type={experiment_type}"
    
    print(f"\n🔍 Fetching ALL available samples for '{experiment_options[experiment_type]}'... (this may take a while)")
    log_message(f"Fetching samples for {experiment_type} started.")

    start_time = time.time()
    processed_count = 0

    while next_page:
        response = requests.get(next_page)
        if response.status_code != 200:
            error_message = f"❌ Error fetching data: {response.status_code}"
            print(error_message)
            log_message(error_message)
            return []
        data = response.json()
        for sample in data["data"]:
            sample_accessions.append(sample["id"])  # 'id' contains the sample accession
            processed_count += 1
            # Log every 1000 samples instead of printing
            if processed_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                log_message(f"Processed {processed_count} samples... (Elapsed Time: {elapsed_time:.2f} seconds)")
        next_page = data["links"].get("next")  # Get next page if available

    total_time = time.time() - start_time
    print(f"\n✅ Retrieved {len(sample_accessions)} samples in {total_time:.2f} seconds.")
    log_message(f"✅ Completed fetching {len(sample_accessions)} samples in {total_time:.2f} seconds.")
    
    return sample_accessions

def fetch_metadata_concurrent(accessions, max_workers=10):
    """Fetch metadata concurrently for each sample and log progress."""
    metadata_list = []
    metadata_url_base = "https://www.ebi.ac.uk/metagenomics/api/latest/samples/"
    
    print("\n🔄 Fetching metadata for each sample concurrently...")
    log_message("Fetching metadata concurrently started.")

    start_time = time.time()

    def get_sample_metadata(accession):
        """Fetch metadata for a single sample with retry for rate limiting."""
        retries = 3
        delay = 5  # seconds to wait if rate limited
        for attempt in range(retries):
            response = requests.get(f"{metadata_url_base}{accession}")
            if response.status_code == 200:
                data = response.json()
                metadata = {
                    "sample_accession": accession,
                    "biome": data["attributes"].get("biome", "N/A"),
                    "environment": data["attributes"].get("environment_material", "N/A"),
                    "temperature": data["attributes"].get("environment_temperature", "N/A"),
                    "salinity": data["attributes"].get("environment_salinity", "N/A"),
                    "pH": data["attributes"].get("environment_ph", "N/A"),
                    "latitude": data["attributes"].get("latitude", "N/A"),
                    "longitude": data["attributes"].get("longitude", "N/A"),
                    "collection_date": data["attributes"].get("collection_date", "N/A"),
                    "study_accession": data["relationships"]["study"]["data"]["id"],
                    "experiment_type": data["attributes"].get("experiment_type", experiment_type)
                }
                return metadata
            elif response.status_code == 429:
                log_message(f"Rate limited for {accession} on attempt {attempt + 1}. Waiting {delay} seconds before retrying.")
                time.sleep(delay)
            else:
                log_message(f"Error {response.status_code} fetching metadata for {accession}.")
                break
        return None

    processed_count = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_accession = {executor.submit(get_sample_metadata, acc): acc for acc in accessions}
        for future in concurrent.futures.as_completed(future_to_accession):
            result = future.result()
            if result:
                metadata_list.append(result)
            processed_count += 1
            if processed_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                log_message(f"Processed metadata for {processed_count} samples... (Elapsed Time: {elapsed_time:.2f} seconds)")

    total_time = time.time() - start_time
    print(f"\n✅ Retrieved metadata for {len(metadata_list)} samples in {total_time:.2f} seconds.")
    log_message(f"✅ Completed fetching metadata for {len(metadata_list)} samples in {total_time:.2f} seconds.")
    
    return metadata_list

# Fetch ALL samples
sample_accessions = fetch_samples(experiment_type)

# Stop if no samples were found
if not sample_accessions:
    print("❌ No samples found for this experiment type.")
    log_message("❌ No samples found. Script terminated.")
    sys.exit(1)

# Fetch metadata concurrently for all sample accessions
metadata = fetch_metadata_concurrent(sample_accessions, max_workers=10)

# Convert metadata to DataFrame
df = pd.DataFrame(metadata)

# Save to CSV with experiment type in filename
csv_filename = f"mgnify_samples_metadata_{experiment_type.replace('-', '_')}.csv"
df.to_csv(csv_filename, index=False)

print(f"\n✅ **Metadata saved to '{csv_filename}'.** 🎉")
log_message(f"✅ Metadata saved to '{csv_filename}'.")

print("\n🎯 Fetching completed. Check `fetch_log.txt` for details.")
log_message("🎯 Fetching completed successfully.")
