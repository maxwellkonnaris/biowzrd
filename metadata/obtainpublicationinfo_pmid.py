import csv
import time
import argparse
from datetime import datetime
import sys
import os
import re

try:
    from Bio import Entrez
except ImportError:
    sys.exit("The Biopython package is required. Install it using `pip install biopython`.")

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
            # Default to comma if detection fails
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
    while True:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            return records
        except Exception as e:
            print(f"Error fetching PMID {pmid}: {e}. Retrying in 30 seconds...")
            time.sleep(30)

def detect_measurement_type(keywords, mesh_terms):
    # Convert all to lowercase for case-insensitive matching
    keywords_str = "; ".join(keywords).lower()
    mesh_str = "; ".join(mesh_terms).lower()
    
    # Check for flow cytometry terms
    flow_terms = ["flow cytometry", "flowcytometry", "flow-cytometry", 
                 "facs", "fluorescence-activated cell sorting"]
    if any(re.search(r'\b' + re.escape(term) + r'\b', keywords_str + mesh_str) for term in flow_terms):
        return "Flow Cytometry"
    
    # Check for PCR terms
    pcr_terms = ["qpcr", "quantitative pcr", "real-time pcr", "rt-pcr", "rt pcr",
                "ddpcr", "digital pcr", "droplet digital pcr",
                "hampcr", "ham pcr"]
    if any(re.search(r'\b' + re.escape(term) + r'\b', keywords_str + mesh_str) for term in pcr_terms):
        for term in pcr_terms:
            if re.search(r'\b' + re.escape(term) + r'\b', keywords_str + mesh_str):
                if term in ["qpcr", "quantitative pcr", "real-time pcr", "rt-pcr", "rt pcr"]:
                    return "qPCR"
                elif term in ["ddpcr", "digital pcr", "droplet digital pcr"]:
                    return "ddPCR"
                elif term in ["hampcr", "ham pcr"]:
                    return "hamPCR"
    
    # Check for spike-in
    if re.search(r'\bspike[\s-]?in\b', keywords_str + mesh_str):
        return "spike-in"
    
    return ""

def detect_sequencing_type(keywords, mesh_terms):
    # Convert all to lowercase for case-insensitive matching
    keywords_str = "; ".join(keywords).lower()
    mesh_str = "; ".join(mesh_terms).lower()
    
    # Check for 16S or amplicon sequencing
    if (re.search(r'\b16s\b', keywords_str + mesh_str) or 
        re.search(r'\bamplicon\b', keywords_str + mesh_str)):
        return "16S/amplicon"
    
    # Check for shotgun sequencing
    if re.search(r'\bshotgun\b', keywords_str + mesh_str):
        return "shotgun"
    
    # Check for metagenomics
    if re.search(r'\bmetagenomic\b', keywords_str + mesh_str):
        return "metagenomics"
    
    return ""

def main():
    args = prompt_if_missing(parse_arguments())

    Entrez.email = args.email
    Entrez.api_key = args.api_key

    pmid_study_pairs = load_pmid_study_pairs(args.pmid_study)

    with open("study_data.csv", mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
            "PMID", "Study Identifier", "Title", "Authors", "Consortium", "Publication Year", "Link",
            "Measurement Type", "Sequencing", "Abstract", "Journal", "Volume", "Issue", "Pages", "DOI", "PMCID",
            "Publication Types", "MeSH Terms", "Keywords", "Affiliations", "Grants", "Retrieved Date"
        ])

        for pmid, study_identifier in pmid_study_pairs:
            today = datetime.today().strftime("%Y-%m-%d")

            if pmid.lower() == "na":
                writer.writerow([pmid, study_identifier] + [""] * 19 + [today])
                continue

            try:
                records = fetch_pubmed_data(pmid)
                article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
                citation = records["PubmedArticle"][0]["MedlineCitation"]
                pubmed_data = records["PubmedArticle"][0].get("PubmedData", {})

                pub_date = article["Journal"]["JournalIssue"]["PubDate"]
                pub_year = pub_date.get("Year", "Unknown")
                title = article.get("ArticleTitle", "")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                abstract_text = article.get("Abstract", {}).get("AbstractText", "")
                if isinstance(abstract_text, list):
                    abstract_text = " ".join(str(a) for a in abstract_text)

                authors = article.get("AuthorList", [])
                author_names = [f"{a['ForeName']} {a['LastName']}" for a in authors if "ForeName" in a and "LastName" in a]

                if study_identifier.upper() == "NA" and authors:
                    first_author_lastname = authors[0].get("LastName", "Unknown")
                    study_identifier = f"{first_author_lastname}{pub_year}"

                journal = article["Journal"].get("Title", "")
                volume = article["Journal"]["JournalIssue"].get("Volume", "")
                issue = article["Journal"]["JournalIssue"].get("Issue", "")
                pages = article.get("Pagination", {}).get("MedlinePgn", "")

                doi = pmcid = ""
                for id_elem in pubmed_data.get("ArticleIdList", []):
                    id_type = id_elem.attributes.get("IdType", "")
                    if id_type == "doi":
                        doi = str(id_elem)
                    elif id_type == "pmc":
                        pmcid = str(id_elem)

                pub_types = [pt for pt in article.get("PublicationTypeList", [])]
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

                consortium = ""
                
                # Detect measurement type and sequencing type
                measurement_type = detect_measurement_type(keywords, mesh_terms)
                sequencing_type = detect_sequencing_type(keywords, mesh_terms)

                writer.writerow([
                    pmid,
                    study_identifier,
                    title,
                    ", ".join(author_names),
                    consortium,
                    pub_year,
                    link,
                    measurement_type,
                    sequencing_type,
                    abstract_text,
                    journal,
                    volume,
                    issue,
                    pages,
                    doi,
                    pmcid,
                    "; ".join(pub_types),
                    "; ".join(mesh_terms),
                    "; ".join(keywords),
                    "; ".join(affiliations),
                    "; ".join(grants),
                    today
                ])

            except Exception as e:
                print(f"Final error for PMID {pmid}: {e}")
                writer.writerow([pmid, study_identifier, "Error retrieving data"] + ["Error"] * 19 + [today])

    print("CSV file 'study_data.csv' has been created.")

if __name__ == "__main__":
    main()
