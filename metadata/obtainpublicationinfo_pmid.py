import csv
import time
import argparse
from datetime import datetime
import sys
import os
import re
from collections import OrderedDict

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
            time.sleep(30)

def detect_measurement_type(keywords, mesh_terms):
    keywords_str = "; ".join(keywords).lower()
    mesh_str = "; ".join(mesh_terms).lower()
    
    flow_terms = ["flow cytometry", "flowcytometry", "flow-cytometry", 
                 "facs", "fluorescence-activated cell sorting"]
    if any(re.search(r'\b' + re.escape(term) + r'\b', keywords_str + mesh_str) for term in flow_terms):
        return "Flow Cytometry"
    
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
    
    if re.search(r'\bspike[\s-]?in\b', keywords_str + mesh_str):
        return "spike-in"
    
    return ""

def detect_sequencing_type(keywords, mesh_terms):
    keywords_str = "; ".join(keywords).lower()
    mesh_str = "; ".join(mesh_terms).lower()
    
    if (re.search(r'\b16s\b', keywords_str + mesh_str) or 
        re.search(r'\bamplicon\b', keywords_str + mesh_str)):
        return "16S/amplicon"
    
    if re.search(r'\bshotgun\b', keywords_str + mesh_str):
        return "shotgun"
    
    if re.search(r'\bmetagenomic\b', keywords_str + mesh_str):
        return "metagenomics"
    
    return ""
    
def get_processed_pmids(output_file):
    processed = set()
    if not os.path.exists(output_file):
        return processed
        
    with open(output_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row and 'PMID' in row:
                processed.add(row['PMID'])
    return processed
    
def get_all_possible_fields(pmid_study_pairs, existing_columns):
    all_fields = set(existing_columns)
    sample_pmids = [p[0] for p in pmid_study_pairs[:3] if p[0].lower() != "na"]
    
    for pmid in sample_pmids:
        try:
            records = fetch_pubmed_data(pmid)
            article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
            citation = records["PubmedArticle"][0]["MedlineCitation"]
            pubmed_data = records["PubmedArticle"][0].get("PubmedData", {})
            
            def add_nested_fields(prefix, data):
                for field, value in data.items():
                    full_field = f"{prefix}_{field}"
                    all_fields.add(full_field)
                    if isinstance(value, dict):
                        add_nested_fields(full_field, value)
                    elif isinstance(value, list) and value and isinstance(value[0], dict):
                        for item in value:
                            add_nested_fields(full_field, item)
            
            add_nested_fields("Article", article)
            add_nested_fields("Citation", citation)
            add_nested_fields("PubmedData", pubmed_data)
                
        except Exception as e:
            print(f"Error processing sample PMID {pmid}: {e}")
            continue
    
    return sorted(all_fields)

def write_with_retry(writer, row_data, output_file, all_columns):
    max_retries = 3
    retry_count = 0
    
    while retry_count < max_retries:
        try:
            # Ensure all fields are present in the row
            for field in all_columns:
                if field not in row_data:
                    row_data[field] = ""
            
            writer.writerow(row_data)
            return all_columns
            
        except ValueError as e:
            if "dict contains fields not in fieldnames" in str(e):
                missing_fields = eval(str(e).split(": ")[1])
                print(f"Adding missing fields to CSV: {missing_fields}")
                
                # Update the columns list
                all_columns.extend(missing_fields)
                all_columns = sorted(set(all_columns))  # Remove duplicates
                
                # Close the current file
                writer.file.close()
                
                # Read existing data
                existing_data = []
                if os.path.exists(output_file):
                    with open(output_file, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f)
                        existing_data = list(reader)
                
                # Rewrite the entire file with new headers
                with open(output_file, 'w', newline='', encoding='utf-8') as f:
                    new_writer = csv.DictWriter(f, fieldnames=all_columns)
                    new_writer.writeheader()
                    for row in existing_data:
                        # Ensure all fields are present in existing rows
                        for field in all_columns:
                            if field not in row:
                                row[field] = ""
                        new_writer.writerow(row)
                
                # Reopen the file in append mode
                file = open(output_file, 'a', newline='', encoding='utf-8')
                writer = csv.DictWriter(file, fieldnames=all_columns)
                
                retry_count += 1
                continue
            else:
                raise
                
    raise ValueError(f"Failed to write after {max_retries} retries")

def main():
    args = prompt_if_missing(parse_arguments())
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    pmid_study_pairs = load_pmid_study_pairs(args.pmid_study)
    output_file = "publication_data.csv"
    
    existing_columns = [
        "PMID", "Study Identifier", "Title", "Authors", "Consortium", "Publication Year", "Link",
        "Measurement Type", "Sequencing", "Abstract", "Journal", "Volume", "Issue", "Pages", "DOI", "PMCID",
        "Publication Types", "MeSH Terms", "Keywords", "Affiliations", "Grants", "Retrieved Date"
    ]
    
    all_columns = get_all_possible_fields(pmid_study_pairs, existing_columns)
    
    if os.path.exists(output_file):
        with open(output_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            existing_headers = next(reader, None)
            if existing_headers:
                all_columns = sorted(set(all_columns) | set(existing_headers))
    
    processed_pmids = get_processed_pmids(output_file)
    file_exists = os.path.exists(output_file)
    
    # Open the file in append mode
    with open(output_file, mode='a' if file_exists else 'w', newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=all_columns)
        if not file_exists:
            writer.writeheader()
            
        new_pmids_processed = 0
        for pmid, study_identifier in pmid_study_pairs:
            if pmid in processed_pmids:
                continue
                
            new_pmids_processed += 1
            row_data = {col: "" for col in all_columns}
            row_data["Study Identifier"] = study_identifier
            row_data["Retrieved Date"] = datetime.today().strftime("%Y-%m-%d")

            if pmid.lower() == "na":
                row_data["PMID"] = pmid
                all_columns = write_with_retry(writer, row_data, output_file, all_columns)
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

                consortium = ""
                
                measurement_type = detect_measurement_type(keywords, mesh_terms)
                sequencing_type = detect_sequencing_type(keywords, mesh_terms)

                row_data.update({
                    "PMID": pmid,
                    "Title": title,
                    "Authors": ", ".join(author_names),
                    "Consortium": consortium,
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
                    "Grants": "; ".join(grants)
                })

                # Add all article fields
                for field in article.keys():
                    field_name = f"Article_{field}"
                    if field_name not in all_columns:
                        all_columns.append(field_name)
                    row_data[field_name] = str(article[field]) if isinstance(article[field], (dict, list)) else article[field]
                
                # Add all citation fields
                for field in citation.keys():
                    field_name = f"Citation_{field}"
                    if field_name not in all_columns:
                        all_columns.append(field_name)
                    row_data[field_name] = str(citation[field]) if isinstance(citation[field], (dict, list)) else citation[field]
                
                # Add all pubmed_data fields
                for field in pubmed_data.keys():
                    field_name = f"PubmedData_{field}"
                    if field_name not in all_columns:
                        all_columns.append(field_name)
                    row_data[field_name] = str(pubmed_data[field]) if isinstance(pubmed_data[field], (dict, list)) else pubmed_data[field]

                all_columns = write_with_retry(writer, row_data, output_file, all_columns)

            except Exception as e:
                print(f"Final error for PMID {pmid}: {e}")
                row_data["PMID"] = pmid
                row_data["Title"] = "Error retrieving data"
                all_columns = write_with_retry(writer, row_data, output_file, all_columns)

    total_pmids = len(pmid_study_pairs)
    skipped_pmids = total_pmids - new_pmids_processed
    print(f"Processing complete. Total PMIDs: {total_pmids}, New: {new_pmids_processed}, Skipped: {skipped_pmids}")
    print(f"Data saved to: {output_file}")

if __name__ == "__main__":
    main()
