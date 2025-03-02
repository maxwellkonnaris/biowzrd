#!/usr/bin/env bash

# Check if the first argument is "--default"
if [ "$1" == "--default" ]; then
  echo "Running in default mode (non-interactive)"
  sequence_type="Amplicon"
  target_gene="16S"
  seq_meth="ion torrent,illumina"
  match="any"
  auto_continue="y"
else
  # Prompt user for options (press Enter to use the default value)
  read -p "Enter sequence type [Amplicon]: " sequence_type
  sequence_type=${sequence_type:-Amplicon}

  read -p "Enter target gene [16S]: " target_gene
  target_gene=${target_gene:-16S}

  read -p "Enter sequencing method (comma separated) [ion torrent,illumina]: " seq_meth
  seq_meth=${seq_meth:-"ion torrent,illumina"}

  read -p "Enter match type (any/all) [any]: " match
  match=${match:-any}

  read -p "When forced retry fails, continue with partial data? [y/N]: " auto_continue
  auto_continue=${auto_continue:-n}
fi

# Desired page size
limit=1000

# URL encode spaces in seq_meth
seq_meth=$(echo "$seq_meth" | sed 's/ /%20/g')

# Build the initial search URL using the provided parameters
page_url="https://api.mg-rast.org/search?sequence_type=${sequence_type}&target_gene=${target_gene}&seq_meth=${seq_meth}&match=${match}&limit=${limit}"

echo "Using search URL: $page_url"

# Download the first chunk
echo "Fetching first page of results from: $page_url"
curl -ks "$page_url" -o "chunk-0.json"

# Parse the JSON
first_json=$(cat "chunk-0.json")

# Get total_count (for reference)
total_count=$(echo "$first_json" | jq '.total_count')
echo "Total results reported by MG-RAST: $total_count"

# Extract 'next' from the first chunk and URL-encode spaces
next_url=$(echo "$first_json" | jq -r '.next // ""' | sed 's/ /%20/g')

# Keep track of how many total records we've downloaded so far
accumulated=0

# Check how many records are in this first chunk
chunk_count=$(echo "$first_json" | jq '.data | length')
accumulated=$((accumulated + chunk_count))

echo "chunk-0.json has $chunk_count records, accumulated total: $accumulated"

# Start the page counter
n=0

# Function to download a page, with up to 3 retries if it's "too small".
# It saves the attempt with the highest record count.
download_chunk_with_retry() {
  local url="$1"
  local outfile="$2"
  local -i attempts=0
  local -i max_attempts=3
  local best_count=0

  while [ $attempts -lt $max_attempts ]; do
    local temp_file
    temp_file=$(mktemp)
    echo "Downloading: $url -> $temp_file (attempt $((attempts+1))/$max_attempts)"
    curl -ks "$url" -o "$temp_file"
    local count
    count=$(jq '.data | length' "$temp_file")
    echo "Chunk has $count records"
    if [ "$count" -gt "$best_count" ]; then
      best_count=$count
      cp "$temp_file" "$outfile"
    fi
    # If chunk has at least 500 records (half of our limit) or zero (last page), we assume it's acceptable.
    if [ "$count" -ge 500 ] || [ "$count" -eq 0 ]; then
      rm "$temp_file"
      return 0
    fi
    rm "$temp_file"
    ((attempts++))
    sleep 30
  done

  echo "⚠️ Gave up after $max_attempts attempts. Best attempt had $best_count records, proceeding with that."
  return 0
}

# While we still have a 'next' URL, keep paginating
while [ -n "$next_url" ]; do
  ((n++))
  outfile="chunk-${n}.json"

  # Download the chunk with retries if it's too small
  download_chunk_with_retry "$next_url" "$outfile"

  # Count how many records we got from the best attempt
  chunk_count=$(jq '.data | length' "$outfile")
  accumulated=$((accumulated + chunk_count))

  echo "chunk-$n.json has $chunk_count records, accumulated total: $accumulated"

  # Extract the new next URL and URL-encode spaces
  new_next=$(jq -r '.next // ""' "$outfile" | sed 's/ /%20/g')

  # If new_next is empty BUT we haven't reached total_count, try one more forced retry.
  if [ -z "$new_next" ] && [ "$accumulated" -lt "$total_count" ]; then
    echo "⚠️ Warning: next is null but we only have $accumulated/$total_count records."
    echo "Trying one more forced retry in case MG-RAST ended early..."
    sleep 30
    download_chunk_with_retry "$next_url" "$outfile"
    chunk_count=$(jq '.data | length' "$outfile")
    accumulated=$((accumulated + chunk_count))
    echo "Forced retry gave us chunk-$n with $chunk_count records, new accumulated total: $accumulated"
    new_next=$(jq -r '.next // ""' "$outfile" | sed 's/ /%20/g')
    if [ -z "$new_next" ]; then
      echo "Still no 'next' URL after forced retry."
      if [ "$auto_continue" = "y" ] || [ "$auto_continue" = "Y" ]; then
        echo "Continuing with partial data as per default setting."
        break
      else
        read -p "Do you want to continue with the partial data (Y) or abort (N)? [Y/n]: " user_choice
        user_choice=${user_choice:-Y}
        if [[ "$user_choice" =~ ^[Yy] ]]; then
          echo "Continuing with partial data."
          break
        else
          echo "Aborting as per user input."
          exit 1
        fi
      fi
    fi
  fi

  # Update next_url for the loop
  next_url="$new_next"
done

# Merge all data into one JSON file
echo "Combining all chunks' data into all_data.json..."
jq -s '[.[].data[]]' chunk-*.json > all_data.json

rm -rf *chunk*
# Final record count
final_count=$(jq 'length' all_data.json)

echo "✅ Download complete!"
echo "We retrieved $final_count records in total (MG-RAST originally reported $total_count)."

jq -r '.[].metagenome_id' all_data.json >> mgrastsamples.txt
echo "Saved sequence IDs to mgrastsamples.txt"
