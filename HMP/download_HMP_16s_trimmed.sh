#!/bin/bash

# Job submission settings (modify as needed for your system)
#SBATCH --job-name=hmp_download
#SBATCH --output=hmp_download.log  # Log file
#SBATCH --error=hmp_download.err   # Error log file
#SBATCH --time=48:00:00            # Max runtime (adjust as needed)
#SBATCH --mem=16G                   # Memory (adjust based on need)
#SBATCH --cpus-per-task=4          # CPU cores required

echo "Starting HMP data download job..."

mkdir -p hmp_16s_trimmed/

aws s3 cp --recursive s3://hmpdcc/hmp1/hhs/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/ihmp/ibdmdb/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/ihmp/momspi/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request 
aws s3 cp --recursive s3://hmpdcc/ihmp/t2d/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request 
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46307/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46309/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46315/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46317/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46319/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46321/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46323/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46327/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46331/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46333/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46335/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46337/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46339/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46877/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA46879/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request
aws s3 cp --recursive s3://hmpdcc/hmp1/demo/PRJNA50637/microbiome/16s/analysis/trimmed/ hmp_16s_trimmed/ --no-sign-request


echo "Download completed!"

