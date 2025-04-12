import pandas as pd
import sys
import argparse

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate input.csv for SLURM pipeline by combining run_accessions.csv and silverman_datarepo.csv.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-r", "--run",
        default="run_accessions.csv",
        help="Path to run_accessions.csv containing Project Accession and Run Accession columns."
    )
    parser.add_argument(
        "-s", "--silverman",
        default="silverman_datarepo.csv",
        help="Path to silverman_datarepo.csv containing accession and sequencingtype columns."
    )
    parser.add_argument(
        "-o", "--output",
        default="input.csv",
        help="Path to output input.csv with Bioproject,RunAccession,SequencingType columns."
    )
    return parser.parse_args()

def main():
    # Parse arguments
    args = parse_arguments()
    run_accessions_file = args.run
    silverman_datarepo_file = args.silverman
    output_file = args.output

    # Read input files
    try:
        run_df = pd.read_csv(run_accessions_file)
        silverman_df = pd.read_csv(silverman_datarepo_file)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: One of the input files is empty")
        sys.exit(1)

    # Validate required columns
    required_run_cols = ["Project Accession", "Run Accession"]
    required_silverman_cols = ["accession", "sequencingtype"]
    if not all(col in run_df.columns for col in required_run_cols):
        print(f"Error: {run_accessions_file} missing required columns: {required_run_cols}")
        sys.exit(1)
    if not all(col in silverman_df.columns for col in required_silverman_cols):
        print(f"Error: {silverman_datarepo_file} missing required columns: {required_silverman_cols}")
        sys.exit(1)

    # Check for duplicate Project Accessions
    duplicate_projects = run_df[run_df["Project Accession"].duplicated(keep=False)]["Project Accession"].unique()
    if len(duplicate_projects) > 0:
        print(f"Warning: Found {len(duplicate_projects)} Project Accessions with multiple Run Accessions: {duplicate_projects[:10]}{'...' if len(duplicate_projects) > 10 else ''}")

    # Get set of valid Project Accessions
    run_accessions = set(run_df["Project Accession"].unique())

    # Filter silverman_datarepo.csv rows
    def has_valid_accession(accession_str):
        if pd.isna(accession_str):
            return False
        # Replace commas with spaces and split
        accessions = accession_str.replace(",", " ").split()
        # Check if any accession matches run_accessions
        return any(acc.strip() in run_accessions for acc in accessions)

    # Apply filter to keep only relevant rows
    silverman_df = silverman_df[silverman_df["accession"].apply(has_valid_accession)]
    if silverman_df.empty:
        print(f"Error: No rows in {silverman_datarepo_file} match Project Accessions from {run_accessions_file}")
        sys.exit(1)

    # Map SequencingType based on case-insensitive substring
    def map_sequencing_type(value):
        if pd.isna(value):
            return None
        value_lower = str(value).lower()
        if "16s" in value_lower:
            return "16S"
        if "shotgun" in value_lower:
            return "meta"
        return None

    silverman_df["sequencingtype"] = silverman_df["sequencingtype"].apply(map_sequencing_type)
    if silverman_df["sequencingtype"].isna().any():
        invalid_types = silverman_df[silverman_df["sequencingtype"].isna()]["accession"].unique()
        print(f"Warning: Skipping entries with unrecognized SequencingType for accession lists: {invalid_types}")
        silverman_df = silverman_df.dropna(subset=["sequencingtype"])

    # Create a mapping of Project Accession to SequencingType
    accession_to_seqtype = {}
    for _, row in silverman_df.iterrows():
        if pd.isna(row["accession"]):
            continue
        # Split accessions
        accessions = row["accession"].replace(",", " ").split()
        seq_type = row["sequencingtype"]
        # Assign sequencingtype to each accession (first match wins)
        for acc in accessions:
            acc = acc.strip()
            if acc in run_accessions and acc not in accession_to_seqtype:
                accession_to_seqtype[acc] = seq_type

    # Add SequencingType to run_df
    run_df["SequencingType"] = run_df["Project Accession"].map(accession_to_seqtype)

    # Check for unmatched Project Accessions
    unmatched = run_df[run_df["SequencingType"].isna()]["Project Accession"].unique()
    if len(unmatched) > 0:
        print(f"Warning: No SequencingType found for Project Accessions: {unmatched}")
        run_df = run_df.dropna(subset=["SequencingType"])

    # Prepare output with correct headers
    output_df = run_df[["Project Accession", "Run Accession", "SequencingType"]].rename(
        columns={
            "Project Accession": "Bioproject",
            "Run Accession": "RunAccession"
        }
    )
    if output_df.empty:
        print(f"Error: No valid data to write to {output_file}")
        sys.exit(1)

    # Verify row count
    input_row_count = len(output_df)
    run_accessions_row_count = len(pd.read_csv(run_accessions_file))
    if input_row_count > run_accessions_row_count:
        print(f"Error: Output has more rows ({input_row_count}) than input ({run_accessions_row_count})")
        sys.exit(1)
    print(f"Input rows: {run_accessions_row_count}, Output rows: {input_row_count}")

    # Write to input.csv
    output_df.to_csv(output_file, index=False)
    print(f"Successfully created {output_file} with {input_row_count} entries")

if __name__ == "__main__":
    main()
