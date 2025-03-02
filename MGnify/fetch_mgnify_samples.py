import requests
import pandas as pd
import argparse
import sys
import time

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
    description="Fetch bulk metadata from MGnify and save it as a CSV file."
)
parser.add_argument(
    "--experiment",
    type=str,
    choices=experiment_options.keys(),
    default="16s-rrna-gene-amplicon",
    help="Specify the experiment type. Available options: metagenomic, 16s-rrna-gene-amplicon (default), 18s-rrna-gene-amplicon, its-gene-amplicon."
)
args = parser.parse_args()
experiment_type = args.experiment

# Clear previous logs
with open(LOG_FILE, "w") as log_file:
    log_file.write("=== Bulk Fetch Script Execution Started ===\n")

print("\n🔹 **MGnify Bulk Metadata Fetcher** 🔹")
print(f"✅ Fetching bulk metadata for: {experiment_options[experiment_type]}")
print(f"📌 Progress is logged to `{LOG_FILE}`.")

def fetch_bulk_metadata(experiment_type):
    """
    Fetch metadata for all samples using the bulk endpoint.
    We assume the API supports a page_size parameter that returns
    all (or most) samples with full details in one (or a few) requests.
    """
    base_url = "https://www.ebi.ac.uk/metagenomics/api/latest/samples"
    # Adjust the page_size as needed. The browser may be using a similar parameter.
    url = f"{base_url}?experiment-type={experiment_type}&page_size=10000"
    metadata_list = []
    start_time = time.time()
    processed_count = 0

    while url:
        response = requests.get(url)
        if response.status_code != 200:
            error_message = f"❌ Error fetching data: {response.status_code}"
            print(error_message)
            log_message(error_message)
            break

        data = response.json()
        for sample in data.get("data", []):
            metadata = {
                "sample_accession": sample.get("id", "N/A"),
                "biome": sample.get("attributes", {}).get("biome", "N/A"),
                "environment": sample.get("attributes", {}).get("environment_material", "N/A"),
                "temperature": sample.get("attributes", {}).get("environment_temperature", "N/A"),
                "salinity": sample.get("attributes", {}).get("environment_salinity", "N/A"),
                "pH": sample.get("attributes", {}).get("environment_ph", "N/A"),
                "latitude": sample.get("attributes", {}).get("latitude", "N/A"),
                "longitude": sample.get("attributes", {}).get("longitude", "N/A"),
                "collection_date": sample.get("attributes", {}).get("collection_date", "N/A"),
                "study_accession": sample.get("relationships", {}).get("study", {}).get("data", {}).get("id", "N/A"),
                "experiment_type": sample.get("attributes", {}).get("experiment_type", experiment_type)
            }
            metadata_list.append(metadata)
            processed_count += 1
            # Optionally log progress every 1000 samples
            if processed_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                log_message(f"Processed {processed_count} samples... (Elapsed Time: {elapsed_time:.2f} seconds)")

        # Move to the next page if available (if pagination is still in effect)
        url = data.get("links", {}).get("next")
    
    total_time = time.time() - start_time
    print(f"\n✅ Retrieved metadata for {len(metadata_list)} samples in {total_time:.2f} seconds.")
    log_message(f"✅ Completed fetching metadata for {len(metadata_list)} samples in {total_time:.2f} seconds.")
    return metadata_list

# Fetch the bulk metadata
metadata = fetch_bulk_metadata(experiment_type)

# Convert metadata to a DataFrame and save to CSV
df = pd.DataFrame(metadata)
csv_filename = f"mgnify_bulk_metadata_{experiment_type.replace('-', '_')}.csv"
df.to_csv(csv_filename, index=False)

print(f"\n✅ **Metadata saved to '{csv_filename}'.** 🎉")
log_message(f"✅ Metadata saved to '{csv_filename}'.")
print("\n🎯 Bulk fetching completed. Check `fetch_log.txt` for details.")
log_message("🎯 Bulk fetching completed successfully.")
