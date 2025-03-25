import os
import argparse

def extract_accession(filename):
    # Extract accession by removing compression/FASTQ suffix
    # and taking the base name (e.g., SRR123456.fastq.gz -> SRR123456)
    if filename.endswith(('.fastq.gz', '.fq.gz')):
        return filename.split('.')[0]
    return None

def create_accession_set(directory):
    accessions = set()
    for filename in os.listdir(directory):
        accession = extract_accession(filename)
        if accession:
            accessions.add(accession)
    return accessions

def save_checkpoint(accessions, output_file):
    with open(output_file, 'w') as f:
        for acc in sorted(accessions):
            f.write(f"{acc}\n")

def main():
    parser = argparse.ArgumentParser(description="Create checkpoint file of accessions in a directory.")
    parser.add_argument("directory", help="Path to directory containing gzipped FASTQ files")
    parser.add_argument("-o", "--output", default="completed_accessions.txt", help="Checkpoint output file (default: completed_accessions.txt)")
    
    args = parser.parse_args()
    
    accessions = create_accession_set(args.directory)
    save_checkpoint(accessions, args.output)
    
    print(f"Checkpoint file saved to {args.output} with {len(accessions)} accessions.")

if __name__ == "__main__":
    main()
