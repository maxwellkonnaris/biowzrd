import csv
import time
import argparse
import shutil
from datetime import datetime
import sys
import os
import re
import tempfile
from collections import OrderedDict

try:
    from Bio import Entrez
except ImportError:
    sys.exit("The Biopython package is required. Install it using `pip install biopython`.")

# Constants
OUTPUT_FILE = "publication_data.csv"

# Keep only the core columns (no dynamic addition)
BASE_COLUMNS = [
    "PMID", "Study Identifier", "Title", "Authors", "Consortium", "Publication Year", "Link",
    "Measurement Type", "Sequencing", "Abstract", "Journal", "Volume", "Issue", "Pages", 
    "DOI", "PMCID", "Publication Types", "MeSH Terms", "Keywords", "Affiliations", 
    "Grants", "Retrieved Date"
]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Fetch PubMed metadata for a list of PMIDs.")
    parser.add_argument("--pmid_study", required=True, help="Path to CSV/TSV/TXT file with PMID, Study Identifier pairs.")
    parser.add_argument("--email", help="NCBI email (required if not set).")
    parser.add_argument("--api_key", help="NCBI API key (optional, but recommended).")
    return parser.parse_args()

def prompt_if_missing(args):
    if not args.email:
        args.email = input("Enter your NCBI email: ")
    if not args.api_key:
        args.api_key = input("Enter your NCBI API key (press Enter to skip): ")
    return args

def detect_delimiter(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        sample = f.read(2048)
        sniffer = csv.Sniffer()
        try:
            dialect = sniffer.sniff(sample)
            return dialect.delimiter
        except csv.Error:
            return ','

def load_pmid_study_pairs(filepath):
    if not os.path.exists(filepath):
        sys.exit(f"File not found: {filepath}")
    
    delimiter = detect_delimiter(filepath)
    pairs = []

    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if len(row) >= 2:
                pairs.append((row[0].strip(), row[1].strip()))
    return pairs

def fetch_pubmed_data(pmid):
    """
    Fetches the data from PubMed for a given PMID.
    Retries every 30 seconds if a 429 (Too Many Requests) or network error occurs.
    """
    while True:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            # Add a small delay to avoid hitting request rate limits too quickly
            time.sleep(0.5)
            return records
        except Exception as e:
            print(f"Error fetching PMID {pmid}: {e}. Retrying in 30 seconds...")
            time.sleep(30)

import re

def detect_measurement_type(keywords, mesh_terms, abstract):
    """
    Returns a recognized measurement type if any known synonyms are found
    in keywords, MeSH terms, or abstract text.
    """
    # Combine all text into one lowercase string
    combined_text = "; ".join(keywords + mesh_terms + [abstract]).lower()
    
    # 1) Flow Cytometry
    flow_synonyms = [
        "flow cytometry", "flowcytometry", "flow-cytometry", "facs",
        "fluorescence-activated cell sorting", "fluorescence activated cell sorting",
        "cytofluorimetric", "cell sorting"
    ]
    if any(re.search(r'\b' + re.escape(term) + r'\b', combined_text) for term in flow_synonyms):
        return "Flow Cytometry"
    
    # 2) qPCR
    qpcr_synonyms = [
        "qpcr", "quantitative pcr", "real-time pcr", "real time pcr", "rt-pcr", "rt pcr",
        "qrt-pcr", "q-rt-pcr", "pcr-based quantification"
    ]
    if any(re.search(r'\b' + re.escape(term) + r'\b', combined_text) for term in qpcr_synonyms):
        return "qPCR"
    
    # 3) ddPCR
    ddpcr_synonyms = [
        "ddpcr", "digital pcr", "droplet digital pcr", "d-pcr"
    ]
    if any(re.search(r'\b' + re.escape(term) + r'\b', combined_text) for term in ddpcr_synonyms):
        return "ddPCR"
    
    # 4) hamPCR
    hampcr_synonyms = [
        "hampcr", "ham pcr"
    ]
    if any(re.search(r'\b' + re.escape(term) + r'\b', combined_text) for term in hampcr_synonyms):
        return "hamPCR"
    
    # 5) Spike-in
    # Looks for "spike-in" or "spike in" (with optional hyphen/whitespace)
    if re.search(r'\bspike[\s-]?in\b', combined_text):
        return "spike-in"
    
    return ""



def detect_sequencing_type(keywords, mesh_terms, abstract):
    """
    Returns a list of recognized sequencing types if any known synonyms are found
    in keywords, MeSH terms, or abstract text.
    """
    combined_text = "; ".join(keywords + mesh_terms + [abstract]).lower()
    detected_types = set()

    type_synonyms = {
        "Amplicon": [
            "amplicon", "amplicon-based", "targeted amplicon", "amplicon sequencing",
            "16s amplicon", "pcr amplicon", "pcr-based", "marker gene sequencing"
        ],
        "Shotgun": [
            "shotgun", "shotgun metagenomic", "whole-genome shotgun", "wgs",
            "shotgun sequencing", "untargeted sequencing"
        ],
        "Metagenomics": [
            "metagenomic", "metagenomics", "whole-metagenome", "metagenomic sequencing",
            "environmental genomics", "microbiome profiling", "meta-genome", "metagenome"
        ],
        "16S": [
            "16s", "16s rrna", "16s rna", "16s ribosomal rna", "16s sequencing",
            "16s gene", "16s-based", "rrs gene"
        ],
        "ITS": [
            "its", "internal transcribed spacer", "its1", "its2", "fungal sequencing",
            "fungal profiling", "its region", "its amplicon"
        ],
        "Metatranscriptomics": [
            "metatranscriptomics", "metatranscriptomic", "rna-based metagenomics",
            "environmental rna", "total rna sequencing", "meta-transcriptome"
        ]
    }

    # Check for all matches
    for seq_type, synonyms in type_synonyms.items():
        if any(re.search(r'\b' + re.escape(term) + r'\b', combined_text) for term in synonyms):
            detected_types.add(seq_type)

    return sorted(detected_types)  


def get_processed_pmids():
    processed = set()
    if not os.path.exists(OUTPUT_FILE):
        return processed
        
    with open(OUTPUT_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row and 'PMID' in row:
                processed.add(row['PMID'])
    return processed

def get_current_headers():
    """
    Returns the headers from the existing CSV if it exists,
    otherwise returns the base columns.
    """
    if not os.path.exists(OUTPUT_FILE):
        return BASE_COLUMNS
        
    with open(OUTPUT_FILE, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        headers = next(reader, None)
        return headers if headers else BASE_COLUMNS

def safe_write(output_file, all_columns, pmid_study_pairs, processed_pmids):
    """
    Writes new rows to a temporary file, then uses shutil.move to replace
    the original CSV. This avoids corruption if the program is interrupted.
    """
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, newline='', encoding='utf-8')
    
    try:
        writer = csv.DictWriter(temp_file, fieldnames=all_columns)
        writer.writeheader()
        
        # If the output file already exists, copy its data first
        if os.path.exists(output_file):
            with open(output_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # Ensure all columns are present
                    for col in all_columns:
                        if col not in row:
                            row[col] = ""
                    writer.writerow(row)
        
        # Process new PMIDs
        new_pmids_processed = 0
        for pmid, study_identifier in pmid_study_pairs:
            if pmid in processed_pmids:
                continue  # Skip already processed PMID
                
            new_pmids_processed += 1
            row_data = {col: "" for col in all_columns}
            row_data["Study Identifier"] = study_identifier
            row_data["Retrieved Date"] = datetime.today().strftime("%Y-%m-%d")

            # Allow "NA" as a special case for PMID
            if pmid.lower() == "na":
                row_data["PMID"] = pmid
                writer.writerow(row_data)
                continue

            try:
                records = fetch_pubmed_data(pmid)
                article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
                citation = records["PubmedArticle"][0]["MedlineCitation"]
                pubmed_data = records["PubmedArticle"][0].get("PubmedData", {})

                # Extract standard fields
                pub_date = article["Journal"]["JournalIssue"]["PubDate"]
                pub_year = pub_date.get("Year", "Unknown")
                title = article.get("ArticleTitle", "")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                abstract_text = article.get("Abstract", {}).get("AbstractText", "")
                if isinstance(abstract_text, list):
                    abstract_text = " ".join(str(a) for a in abstract_text)

                authors = article.get("AuthorList", [])
                author_names = [
                    f"{a['ForeName']} {a['LastName']}"
                    for a in authors if "ForeName" in a and "LastName" in a
                ]

                # If "Study Identifier" is NA, try to set it from first author + year
                if study_identifier.upper() == "NA" and authors:
                    first_author_lastname = authors[0].get("LastName", "Unknown")
                    study_identifier = f"{first_author_lastname} {pub_year}"

                journal = article["Journal"].get("Title", "")
                volume = article["Journal"]["JournalIssue"].get("Volume", "")
                issue = article["Journal"]["JournalIssue"].get("Issue", "")
                pages = article.get("Pagination", {}).get("MedlinePgn", "")

                doi = ""
                pmcid = ""
                for id_elem in pubmed_data.get("ArticleIdList", []):
                    id_type = id_elem.attributes.get("IdType", "")
                    if id_type == "doi":
                        doi = str(id_elem)
                    elif id_type == "pmc":
                        pmcid = str(id_elem)

                pub_types = [str(pt) for pt in article.get("PublicationTypeList", [])]
                mesh_terms = [str(m["DescriptorName"]) for m in citation.get("MeshHeadingList", [])]
                keyword_list = citation.get("KeywordList", [])
                keywords = [str(kw) for group in keyword_list for kw in group]

                affiliations = []
                for author in authors:
                    for aff in author.get("AffiliationInfo", []):
                        affiliations.append(aff.get("Affiliation", ""))

                grants = []
                for grant in citation.get("GrantList", []):
                    agency = grant.get("Agency", "")
                    grant_id = grant.get("GrantID", "")
                    grants.append(f"{agency} ({grant_id})" if grant_id else agency)

                measurement_type = detect_measurement_type(keywords, mesh_terms, abstract_text)
                sequencing_type = detect_sequencing_type(keywords, mesh_terms, abstract_text)

                # Populate row data
                row_data.update({
                    "PMID": pmid,
                    "Title": title,
                    "Authors": ", ".join(author_names),
                    "Publication Year": pub_year,
                    "Link": link,
                    "Measurement Type": measurement_type,
                    "Sequencing": sequencing_type,
                    "Abstract": abstract_text,
                    "Journal": journal,
                    "Volume": volume,
                    "Issue": issue,
                    "Pages": pages,
                    "DOI": doi,
                    "PMCID": pmcid,
                    "Publication Types": "; ".join(pub_types),
                    "MeSH Terms": "; ".join(mesh_terms),
                    "Keywords": "; ".join(keywords),
                    "Affiliations": "; ".join(affiliations),
                    "Grants": "; ".join(grants),
                    "Study Identifier": study_identifier
                })

                writer.writerow(row_data)

            except Exception as e:
                print(f"Error processing PMID {pmid}: {e}")
                row_data["PMID"] = pmid
                row_data["Title"] = f"Error: {str(e)[:100]}"
                writer.writerow(row_data)

        # Close the temp file so contents are flushed
        temp_file.close()

        # Replace the original file using shutil.move (works across devices)
        if os.path.exists(output_file):
            os.remove(output_file)
        shutil.move(temp_file.name, output_file)
        
        return new_pmids_processed

    except Exception as e:
        print(f"Critical error during processing: {e}")
        if os.path.exists(temp_file.name):
            os.remove(temp_file.name)
        raise

def main():
    args = prompt_if_missing(parse_arguments())
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    pmid_study_pairs = load_pmid_study_pairs(args.pmid_study)
    processed_pmids = get_processed_pmids()
    all_columns = get_current_headers()

    batch_size = 100
    total_pmids = len(pmid_study_pairs)
    total_processed = 0
    
    for i in range(0, len(pmid_study_pairs), batch_size):
        batch = pmid_study_pairs[i:i+batch_size]
        new_processed = safe_write(OUTPUT_FILE, all_columns, batch, processed_pmids)
        total_processed += new_processed
        processed_pmids.update(pmid for pmid, _ in batch if pmid.lower() != "na")
        
        # Update headers in case new fields were added (in this version we do NOT add new fields)
        all_columns = get_current_headers()
        
        print(f"Processed batch {i//batch_size + 1}: {new_processed} new records")

    skipped_pmids = total_pmids - total_processed
    print(f"\nProcessing complete. Total PMIDs: {total_pmids}, New: {total_processed}, Skipped: {skipped_pmids}")
    print(f"Data saved to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
