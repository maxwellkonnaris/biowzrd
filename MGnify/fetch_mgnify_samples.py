import requests
import pandas as pd
import argparse
import sys

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

# Add argument for experiment type
parser.add_argument(
    "--experiment",
    type=str,
    choices=experiment_options.keys(),
    default="16s-rrna-gene-amplicon",
    help="Specify the experiment type. Available options: metagenomic, 16s-rrna-gene-amplicon (default), 18s-rrna-gene-amplicon, its-gene-amplicon."
)

# Parse command-line arguments
args = parser.parse_args()
experiment_type = args.experiment  # Get the selected experiment type

# Print usage instructions for clarity
print("\n🔹 **MGnify Metadata Fetcher** 🔹")
print("This script retrieves metadata from MGnify for a selected experiment type and saves it as a CSV file.\n")
print("📌 **Available Experiment Types:**")
for key, value in experiment_options.items():
    print(f"  {key} - {value}")

print("\n✅ **Using Experiment Type:**", experiment_options[experiment_type])
print("\n📌 **Example Usage:**")
print("   python fetch_mgnify_samples.py --experiment metagenomic")
print("   python fetch_mgnify_samples.py --experiment 18s-rrna-gene-amplicon\n")

# Base URL for fetching samples
base_url = "https://www.ebi.ac.uk/metagenomics/api/latest/samples"

# Function to fetch sample accessions
def fetch_samples(experiment_type):
    sample_accessions = []
    next_page = f"{base_url}?experiment-type={experiment_type}"
    
    print(f"\n🔄 Fetching samples for '{experiment_options[experiment_type]}'...")
    
    while next_page:
        response = requests.get(next_page)
        
        # Error handling for bad requests
        if response.status_code != 200:
            print(f"❌ Error fetching data: {response.status_code}")
            sys.exit(1)
        
        data = response.json()

        # Extract sample accessions
        for sample in data["data"]:
            sample_accessions.append(sample["id"])  # 'id' contains the sample accession
        
        # Get next page if available
        next_page = data["links"].get("next")
    
    print(f"✅ Retrieved {len(sample_accessions)} sample accessions.\n")
    return sample_accessions

# Function to fetch metadata for samples
def fetch_metadata(accessions):
    metadata_list = []
    metadata_url_base = "https://www.ebi.ac.uk/metagenomics/api/latest/samples/"
    
    print("🔄 Fetching metadata for each sample...")

    for accession in accessions:
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
            metadata_list.append(metadata)
    
    print(f"✅ Retrieved metadata for {len(metadata_list)} samples.\n")
    return metadata_list

# Fetch samples based on user-selected experiment type
sample_accessions = fetch_samples(experiment_type)

# Stop if no samples were found
if not sample_accessions:
    print("❌ No samples found for this experiment type.")
    sys.exit(1)

metadata = fetch_metadata(sample_accessions)

# Convert metadata to DataFrame
df = pd.DataFrame(metadata)

# Save to CSV with experiment type in filename
csv_filename = f"mgnify_samples_metadata_{experiment_type.replace('-', '_')}.csv"
df.to_csv(csv_filename, index=False)

print(f"✅ **Metadata saved to '{csv_filename}'.** 🎉")
