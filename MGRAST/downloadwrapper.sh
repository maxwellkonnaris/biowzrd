#!/bin/bash
# wrapper.sh

# Calculate the total number of samples
NUM_SAMPLES=$(wc -l < mgrastsamples.txt)
echo "Submitting job array with ${NUM_SAMPLES} tasks"

# Submit the download script with the calculated array range (50 concurrent jobs)
sbatch --array=1-${NUM_SAMPLES}%50 download_mgrast.sh

