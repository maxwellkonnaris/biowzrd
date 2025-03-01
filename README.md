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
