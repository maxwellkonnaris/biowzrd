## Downloading files on high performance compute cluster (downloadingfiles/)
Packages required: NCBI toolkit, fastq-dl
1. create txt file of NCBI or ENA study BioProject accessions called studies.txt
2. run binary file obtainstudyaccessions.sh
   - this will give you a list of all of the run accessions in the bioproject
3. run binary file downloadingfastq.sh

## Quality control (nonassemblyQC/ or assemblyQC/)
Packages required: fastp
1. Perform absolute minimal quality control from publicly available sequences. This is with the assumption that host reads and adapters, if not paired-end sequences, are removed. In NCBI SRA see this statement: https://www.ncbi.nlm.nih.gov/sra/docs/submit/ for Metagenomic data. 

## Kmer based analysis (kmermining/)
Packages required: Jellyfish
1. Perform QC from nonassemblyQC/minimalqualitycontrol.sh. This is done with intention for non assembly based downstream tasks.
2. Specify a range of Kmer lengths to count and mine from each read in each trimmed fastq file after running minimalqualitycontrol.sh 

## Human Microbiome Project (HMP/)
Packages required: aws (see download script for linux)
1. Navigate to home directory and Install aws (installaws.sh)
2. Download all trimmed 16s sequence files in fa.bzip format (downloadHMP_16s_trimmed.sh)
3. Remove reads with length < 30 (removereads.sh)

## MG-RAST (MGRAST/)
1. Download the metadata to obtain metagenome_ids. You can specify the query further as prompted or add the --default flag to download amplicon 16s samples from ion torrent and illumina sequencing technology. (sbatch retrievemetadata.sh)
2. sh downloadmgrast.sh is the combined script which will submit slurm jobs. This will begin downloading preprocessed and host removed fasta files for the metagenome_ids you've specified, in fasta.gz format. 

## MGnify (MGnify/)
1. Download the metadata to obtain sample ids. You can specify the query further to obtain either shotgun metagenomic, 16s/18s/ITS amplicon sequence metadata which is outputted in CSV format. (python fetch_mgnify_samples.py)
2. Use the sample accession to download the sequencing files.
