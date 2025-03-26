## Setup

Follow these steps to **clone the repository** and **run the setup script** to configure your environment.
1. Build conda environment and make sure git-lfs is installed, run:
   ```bash
   conda create -n biowzrd
   conda activate biowzrd
   conda install git-lfs -y
   ```
3. Clone the Repository from GitHub, run:
   ```bash
   git clone git@github.com:maxwellkonnaris/biowzrd.git
   cd biowzrd
   ```
4. Update conda environment for dependencies and required packages, run:
   ```bash
   conda env update --name biowzrd --file environment.yml --yes
   ```
   or, if you are working with micromamba instead
   ```bash
   micromamba update --name biowzrd --file environment.yml --yes
   ```
5. Set permissions for setup.py and Add executable permissions to all script, run:
   ```bash
   chmod +x setup.py
   python setup.py
   ```

## Papers (papers/)

```
📦 papers
├── 📂 papers
│   ├── 📄 {year}_{journal}_{author}_{description}.pdf
├── 📄 .gitattributes
├── 📄 .gitmodules
├── 📄 bibliography.bib - contains BibTeX for all papers with a note field that contains keywords
├── 📄 doi_links.txt - contains DOI links for all papers paired with keywords
├── 📄 run_bib.py - creates and adds BibTeX to bibliography.bib based on doi_links.txt
├── 📄 search_bib.py - functionality to search bibliography.bib based on `{year}, {title}, {keywords}, {author}`
```

Adding to bibliography.bib:

1. Find paper
2. Add DOI and keywords to doi\_links.txt 
    - structure: <doi link> - keyword1, keyword2
```
python run_bib.py
```

Searching bibliography.bib:

```
python search_bib.py "keyword1 keyword2" 
```

## Downloading files on high performance compute cluster (downloadingfiles/)
Packages required: biopython

1. create txt file of NCBI or ENA study BioProject accessions called studies.txt
2. run binary file fetch_runaccession.py
   - this will give you a list of all of the run accessions in the bioproject

```text
             +-----------------------------+
             |       studies.txt           |
             +-----------------------------+
                        |
                        v
             +-----------------------------+
             |    fetch_runaccession.py    |
             +-----------------------------+
                        |
                        v
             +-----------------------------+
             |    run_accessions.txt       |
             +-----------------------------+
                        |
                        v
        +-------------------------------------------------+
        |             execute_slurmmultijob.sh            |
        +-------------------------------------------------+
          |         ^         |                |         ^ 
          |         |         |                |         |
   submit_job  cleanup   manage tokens         |         | _ _ _ _ _ _
         |                                     |                       |
         v                                     |                       |
     +-----------------------------+           |       +-----------------------------+  
     | Job Script for Each ITEM    |        COMPLETE   |    check_accessions.py      | 
     +-----------------------------+           |       +-----------------------------+     
             |                                 |                      ^
             v                                 |                      |
   +-------------------------------+           |        +-----------------------------+ 
   |        download_fastq.py      |           |        |  check_doublecheckfastq.py  |
   +-------------------------------+           |        +-----------------------------+ 
      | env vars (e.g. ACCESSION)              v        
      | download + fallback logic      +-------------------------------+                
      | metadata processing            | Job Script for remaining .SRA |
      | token release (atexit)         +-------------------------------+ 
      | checkpoint logging                
```
3. run binary file execute_slurmmultijob.sh passing in the arguments below as an example:

```bash
./execute_slurmmultijob.sh \
  -i run_accessions.txt \
  -o completed_accessions.txt \
  -m 20 \
  -p fastq_ \
  -C "python downloadingfiles/download_fastq.py" \
  -T "04:00:00" \
  -M "8G" \
  --api-key "your-ncbi-api-key" \
  --email "you@example.com" \
  --export "WORKDIR=$(pwd); CHECKPOINT_FILE=$(pwd)/completed_accessions.txt; CHECKPOINT_LOCK_FILE=$(pwd)/checkpoint.lock; DEBUG_LOCK=$(pwd)/debug.lock; TOKEN_FILE=$(pwd)/.job_tokens; TOKEN_LOCK_FILE=$(pwd)/.job_tokens.lock; COMBINED_METADATA=$(pwd)/combined_metadata.tsv"

```

4. Alternatively, if you want/have access to the raw fastq files. I've created a web scraper that works with the same output from obtainstudyaccessions.sh. First create a conda environment:
```bash
conda create -n sra_downloader -y
conda activate sra_downloader
conda install -y selenium requests beautifulsoup4 tqdm firefox geckodriver sra-tools entrez-direct enaBrowserTools
conda install -y -c conda-forge awscli wget
```
5. run python obtainrawfastqdownloadlinks.py
   - this will give you the accession downloadlink in a tab separated txt file format. 
7. run python download_fastq.py
   - you wont be able to successfully download the raw fastq files unless you have access. some are publicly available, but most are not and require a aws key which can be added to aws configure.
   - this will download the raw fastq files from the links scraped with obtainrawfastqdownloadlinks.py

If you wanted to obtain the full SRA metadata with only runinfo for 16s amplicon sequences from Illumina sequencers:
```bash
esearch -db sra -query "16S[All Fields] AND amplicon[Strategy] and Illumina[Platform]" | efetch -format runinfo > sra_16s_metadata.csv
```

The -format runinfo command should return a table with columns like:

- Run (Run accession ID)
- ReleaseDate (Date when data was made public)
- LoadDate (Date when data was loaded into the SRA database)
- Spots (Number of sequencing reads)
- Bases (Total number of bases)
- Experiment (Experiment accession ID)
- Sample (Sample accession ID)
- BioProject (BioProject ID)
- BioSample (BioSample ID)
- LibraryStrategy (e.g., AMPLICON)
- LibrarySource (e.g., METAGENOMIC)
- LibrarySelection (e.g., PCR)
- LibraryLayout (SINGLE or PAIRED)
- Platform (e.g., ILLUMINA)
- Instrument (e.g., Illumina MiSeq)
- Study (Study accession ID)
- CenterName (Name of the sequencing center)

## File manangement (filemanagement/)
Packages required: GNU parallel
1. Transfer an entire large directory in chunks with paralellel processing - transferdirectory.sh

```bash
# If you would like to track the live progress
watch tail -n 20 transfer_*.log  
```

## Quality control (nonassemblyQC/ or assemblyQC/)
Packages required: fastp, BBMap suite
- nonassemblyQC/
1. Perform absolute minimal quality control from publicly available sequences. This is with the assumption that host reads and adapters, if not paired-end sequences, are removed. This script uses the BBMap adapter reference for commonly used Illumina adapters on the 3' end of a sequence and removes them. Uses BBDuk to mask bp with PHRED < 20 and minimal length of read = 50. - qc_minimal.sh. In NCBI SRA see this statement: https://www.ncbi.nlm.nih.gov/sra/docs/submit/ for Metagenomic data.
2. If using downstream sequence count taxonomic classification like DADA2 or mOTUs or MetaPhlAn. -qc_taxacounts.sh
- assemblyQC/

## Kmer based analysis (kmermining/)
Packages required: Jellyfish
1. Perform QC from nonassemblyQC/minimalqualitycontrol.sh. This is done with intention for non assembly based downstream tasks.
2. Specify a range of Kmer lengths to count and mine from each read in each trimmed fastq file after running minimalqualitycontrol.sh

## Metadata (metadata/)
```bash
./merge_fastq_metadata.py -r output.tsv -s combined_metadata.tsv -o merged_output.tsv
```

## Microbiome Count Tables (microbiomecounts/)
```bash
# Download the rdp classifier, check updates here: https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/

wget https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.14.zip/download -O rdp_classifier_2.14.zip
wget https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/releaseNotes/release_2.14_note.txt/download -O release_2.14_note.txt
unzip rdp_classifier_2.14.zip
cd rdp_classifier_2.14
```
OR to use the training data in dada2:
```bash
wget https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/RDPClassifier_16S_trainsetNo19_QiimeFormat.zip/download -O rdp_trainset.zip
unzip rdp_trainset.zip
```

```bash
## For Amplicon:

# For QC with Fastp and DADA2 classifier (RDP classifier or SILVA classifier):
conda create -dada2 -y -c bioconda -c conda-forge fastp
conda activate dada2
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("dada2", ask = FALSE)'

## For Metagenomics:

# For QC with Fastp and MetaPhlAn4 classifier:
conda create -n metaphlan -y -c bioconda -c conda-forge metaphlan bowtie2 fastp
conda activate metaphlan
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103  

# For QC with Fastp and mOTU classifier:
conda create -n motus -y -c bioconda motus fastp
conda activate motus
motus downloadDB  
```

## List of studies with external total scale measurements paired with sequencing data
Visit: https://docs.google.com/spreadsheets/d/13b4Toscse0MjyAGYt1zfWoPxSRpyuvENHVGKBLYwAAw/edit?usp=sharing

## Human Microbiome Project (HMP/)
Packages required: aws (see download script for linux)
1. To download via aws, navigate to home directory and Install aws (installaws.sh)
2. Download all trimmed 16s sequence files in fa.bzip format (downloadHMP_16s_trimmed.sh)
3. Remove reads with length < 30 (removereads.sh)
4. Metadata or downloading via SRA:

To see how many sequences there are for HMP amplicon sequences generated from illumina sequencers
```R
#library(devtools)
#install_github("ropensci/rentrez")
library(rentrez)
search_results <- entrez_search(db="sra", term="Human Microbiome Project AND amplicon[Strategy] AND Illumina[Platform]", retmax=0, use_history=TRUE)
```

Or instead to jump right to downloading the metadata. Note: Default without the -q (query flag) is what is given below.
```bash
# Download all of the metadata for a specific query:

./querySRA.sh -q "Human Microbiome Project AND 16S[Strategy] AND Illumina[Platform]"
```

## MG-RAST (MGRAST/)
1. Download the metadata to obtain metagenome_ids. You can specify the query further as prompted or add the --default flag to download amplicon 16s samples from ion torrent and illumina sequencing technology. (sbatch retrievemetadata.sh)
2. sh downloadmgrast.sh is the combined script which will submit slurm jobs. This will begin downloading preprocessed and host removed fasta files for the metagenome_ids you've specified, in fasta.gz format. 

## MGnify (MGnify/)
Packages optional: git lfs
1. Download the metadata to obtain sample ids. You can specify the query further to obtain either shotgun metagenomic, 16s/18s/ITS amplicon sequence metadata which is outputted in CSV format. (python fetch_mgnify_samples.py) ----- metadata file provided for 16s amplicon sequences, but need to pull with git lfs pull
2. Obtain the unique accessions or sample IDs to generate the output for downloadingfastq.sh

```bash
awk -F',' 'NR>1 && NF > 1 && !seen[$1]++ && $1 ~ /^[A-Za-z0-9_-]+$/ {print $5}' mgnify_samples_16s-rrna-gene-amplicon.csv > run_accessions.txt
```
   
4. Use the sample accession to download the sequencing files with downloadingfastq.sh
5. If youd like other information about each study, you can use the "study_ids" column. Heres an example where MGYS00006745 represents the study_id. Just copy this into your browser or automate this with a script:

```code
https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00006745 
```

## Plotting (plotting/):

1. Color palletes: https://www.simplifiedsciencepublishing.com/resources/best-color-palettes-for-scientific-figures-and-data-visualizations
2. Color blindness viewer: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
3. In R, the tidyplots package is nice (https://tidyplots.org/), but if making complicated larger figures use ggplot and a custom theme like in theme.R. Ive also listed my 4 color colorway. You can load this into your R with:

```R
source("theme.R")
```

## Garbage, but got it done (garbage/)

This is literally a folder full of garbage scripts that I used to do select analysis. A lot of this code is chatgpt generated (not all) because I just needed to get something done quickly. Was not made for users.
