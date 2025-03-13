import pandas as pd
import argparse
import os

def format_row_as_text(row):
    """
    Convert a Pandas Series (row) into a structured text prompt.
    """
    return " | ".join([f"{col}: {row[col]}" for col in row.index])

def convert_table_to_text(input_file, output_file):
    """
    Reads a CSV file, converts each row to a structured text format, and saves to a file.

    :param input_file: Path to the input CSV file.
    :param output_file: Path to save the formatted text output.
    """
    try:
        # Load dataset
        df = pd.read_csv(input_file)

        # Convert each row to text format
        formatted_rows = df.apply(format_row_as_text, axis=1).tolist()

        # Save to a text file
        with open(output_file, "w", encoding="utf-8") as f:
            for row_text in formatted_rows:
                f.write(row_text + "\n")

        print(f"✅ Formatted text prompts saved to: {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert a tabular CSV dataset into structured text format.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_file", nargs="?", default="formatted_prompts.txt",
                        help="Path to save the formatted text output (default: formatted_prompts.txt)")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"❌ Error: The file '{args.input_file}' does not exist.")
    else:
        convert_table_to_text(args.input_file, args.output_file)
