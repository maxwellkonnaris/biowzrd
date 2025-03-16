#!/usr/bin/env python3

import pandas as pd
import argparse

def merge_files(read_level_file, sample_level_file, output_file, chunksize=100000):
    """
    Merges read-level FASTQ metadata with sample-level metadata from runinfo in SRA.

    Args:
        read_level_file (str): Path to the read-level metadata TSV file.
        sample_level_file (str): Path to the sample-level metadata TSV file.
        output_file (str): Path to save the merged output TSV file.
        chunksize (int): Number of rows to process at a time to handle large files.
    """
    print(f"Loading sample metadata from {sample_level_file}...")
    
    # Load sample metadata into a dictionary for fast lookups
    sample_metadata = pd.read_csv(sample_level_file, sep="\t", index_col="Run").to_dict(orient="index")

    print(f"Processing read-level metadata from {read_level_file} in chunks...")
    
    chunks = []
    
    with pd.read_csv(read_level_file, sep="\t", chunksize=chunksize) as reader:
        for i, chunk in enumerate(reader):
            print(f"Processing chunk {i+1}...")
            
            # Map 'Accession' column in read-level file to 'Run' in sample metadata
            chunk["Sample_Metadata"] = chunk["Accession"].map(sample_metadata)
            
            # Expand dictionary values into separate columns
            sample_metadata_df = pd.DataFrame(chunk["Sample_Metadata"].tolist(), index=chunk.index)
            
            # Merge the extracted metadata with the chunk
            chunk = pd.concat([chunk, sample_metadata_df], axis=1)
            chunk.drop(columns=["Sample_Metadata"], inplace=True)
            
            chunks.append(chunk)
    
    print(f"Writing merged data to {output_file}...")
    merged_df = pd.concat(chunks, ignore_index=True)
    merged_df.to_csv(output_file, sep="\t", index=False)

    print("Merging complete. Output saved to", output_file)

if __name__ == "__main__":
    # Argument parsing for command-line execution
    parser = argparse.ArgumentParser(description="Merge FASTQ read-level metadata with sample metadata efficiently.")
    
    parser.add_argument("-r", "--read", required=True, help="Path to read-level metadata TSV file (output.tsv).")
    parser.add_argument("-s", "--sample", required=True, help="Path to sample-level metadata TSV file (combined_metadata.tsv).")
    parser.add_argument("-o", "--output", required=True, help="Path to save merged output TSV file (merged_output.tsv).")
    
    args = parser.parse_args()
    
    # Call merge function with provided arguments
    merge_files(args.read, args.sample, args.output)
