import pandas as pd
import argparse
import os
import requests
import json
import time
from tqdm import tqdm

# Import local model only if needed (avoids unnecessary dependencies)
try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    SentenceTransformer = None  # Avoids errors if not using local mode

"""
##########################################################################################
#                                      EMBEDDING SCRIPT                                  #
##########################################################################################
# This script allows users to embed tabular data (CSV/TSV) using either:                 #
#   1. The Hugging Face API (default)                                                    #
#   2. A local model (via `sentence-transformers`)                                       #
#                                                                                        #
# ================================ HOW TO USE THE SCRIPT =============================== #
#                                                                                        #
# ▶ Run with Hugging Face API (Default)                                                  #
#   python embed_csv_huggingface.py your_data.csv your_embeddings.csv                    #
#                                                                                        #
# ▶ Use a Local Model Instead                                                            #
#   python embed_csv_huggingface.py your_data.csv your_embeddings.csv --use_local        #
#   (Requires `sentence-transformers`: Install with `pip install sentence-transformers`) #
#                                                                                        #
# ▶ Change the Model                                                                     #
#   python embed_csv_huggingface.py your_data.csv your_embeddings.csv --model sentence-transformers/all-mpnet-base-v2 #
#                                                                                        #
# ▶ Adjust Batch Size for API (Default: 16)                                              #
#   python embed_csv_huggingface.py your_data.csv your_embeddings.csv --batch_size 32    #
#                                                                                        #
##########################################################################################
"""

def detect_file_format(input_file):
    """
    Detect whether the input file is a CSV or TSV based on its extension.
    """
    if input_file.endswith(".tsv"):
        return "tsv"
    elif input_file.endswith(".csv"):
        return "csv"
    else:
        raise ValueError("Unsupported file format. Please use a .csv or .tsv file.")

def format_dataframe_as_text(df):
    """
    Convert a Pandas DataFrame into a list of structured text prompts.
    """
    return df.apply(lambda row: " | ".join([f"{col}: {row[col]}" for col in row.index]), axis=1).tolist()

def get_embeddings_api(texts, model_name, api_key, max_retries=5, initial_delay=2):
    """
    Calls Hugging Face API to get embeddings for a batch of texts with retry logic.
    """
    API_URL = f"https://api-inference.huggingface.co/pipeline/feature-extraction/{model_name}"
    headers = {"Authorization": f"Bearer {api_key}"}
    retries = 0
    delay = initial_delay

    while retries < max_retries:
        try:
            response = requests.post(API_URL, headers=headers, json={"inputs": texts})
            
            if response.status_code == 429:  # Rate limit exceeded
                print(f"⚠️ Rate limit hit! Retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
                retries += 1
                continue

            response.raise_for_status()
            return response.json()

        except requests.exceptions.RequestException as e:
            print(f"❌ API request failed: {e}")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
            retries += 1

    print("❌ Failed to retrieve embeddings after multiple retries.")
    return [None] * len(texts)  # Return None for each failed row to maintain DataFrame integrity

def get_embeddings_local(texts, model_name):
    """
    Generates embeddings using a local model from sentence-transformers.
    """
    if SentenceTransformer is None:
        raise ImportError("sentence-transformers is not installed. Install it using `pip install sentence-transformers`.")

    model = SentenceTransformer(model_name)
    return model.encode(texts).tolist()

def process_file(input_file, output_file, model_name, use_local, api_key, batch_size=16):
    """
    Reads a CSV or TSV file, converts each row to structured text, generates embeddings in batches,
    and saves to a new file while handling API rate limits or using a local model.
    """
    try:
        # Detect file format
        file_format = detect_file_format(input_file)

        # Load dataset with appropriate delimiter
        delimiter = "\t" if file_format == "tsv" else ","
        df = pd.read_csv(input_file, delimiter=delimiter)

        # Convert rows to structured text
        texts = format_dataframe_as_text(df)

        # Process embeddings using API or local model
        embeddings = []
        if use_local:
            print("🔄 Using local model for embeddings...")
            embeddings = get_embeddings_local(texts, model_name)
        else:
            print("🌍 Using Hugging Face API for embeddings...")
            for i in tqdm(range(0, len(texts), batch_size), desc="Embedding batches"):
                batch_texts = texts[i:i+batch_size]
                batch_embeddings = get_embeddings_api(batch_texts, model_name, api_key)

                if batch_embeddings is None:
                    print("⚠️ API call failed, reducing batch size...")
                    batch_size = max(1, batch_size // 2)  # Reduce batch size and retry
                    continue

                embeddings.extend(batch_embeddings)

        # Add embeddings to the DataFrame
        df['embedding'] = embeddings

        # Save results with the correct format
        output_delimiter = "\t" if file_format == "tsv" else ","
        df.to_csv(output_file, index=False, sep=output_delimiter)

        print(f"✅ Embeddings saved to: {output_file}")

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Embed tabular CSV/TSV data using Hugging Face API or a local model.")
    parser.add_argument("input_file", help="Path to the input CSV/TSV file")
    parser.add_argument("output_file", nargs="?", default="embedded_data.csv", help="Path to save the embeddings (default: embedded_data.csv)")
    parser.add_argument("--model", default="sentence-transformers/all-MiniLM-L6-v2", help="Hugging Face model or local model to use for embeddings")
    parser.add_argument("--use_local", action="store_true", help="Use a local model instead of the Hugging Face API")
    parser.add_argument("--api_key", default=os.getenv("HUGGINGFACE_API_KEY"), help="Hugging Face API Key (or set as env variable)")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size for API requests (default: 16)")

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"❌ Error: The file '{args.input_file}' does not exist.")
    elif not args.use_local and not args.api_key:
        print("❌ Error: Hugging Face API key is missing. Pass it via --api_key or set as an environment variable.")
    else:
        process_file(args.input_file, args.output_file, args.model, args.use_local, args.api_key, args.batch_size)
