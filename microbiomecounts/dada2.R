#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript run_dada2.R <forward_fastq> [<reverse_fastq>]", call.=FALSE)
}

# Detect if single-end or paired-end
single_end <- length(args) == 1
input_F <- args[1]
input_R <- if (!single_end) args[2] else NULL

# Output filenames
output_asv <- paste0("dada2_results_", gsub(".fastq.gz", "", basename(input_F)), ".rds")
output_tax <- paste0("dada2_taxonomy_", gsub(".fastq.gz", "", basename(input_F)), ".rds")

# Define paths
filtered_path <- "filtered_fastq/"
if (!dir.exists(filtered_path)) dir.create(filtered_path)

filt_F <- file.path(filtered_path, paste0("filtered_", basename(input_F)))
filt_R <- if (!single_end) file.path(filtered_path, paste0("filtered_", basename(input_R))) else NULL

# Filtering and trimming
if (single_end) {
  out <- filterAndTrim(input_F, filt_F, truncLen=250, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
} else {
  out <- filterAndTrim(c(input_F, input_R), c(filt_F, filt_R), truncLen=c(250,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
}

# Learn error rates
errF <- learnErrors(filt_F, multithread=TRUE)
if (!single_end) {
  errR <- learnErrors(filt_R, multithread=TRUE)
}

# Dereplication
derep_F <- derepFastq(filt_F)
names(derep_F) <- gsub(".fastq.gz", "", basename(filt_F))

if (!single_end) {
  derep_R <- derepFastq(filt_R)
  names(derep_R) <- gsub(".fastq.gz", "", basename(filt_R))
}

# DADA2 core processing
dada_F <- dada(derep_F, err=errF, multithread=TRUE)

if (single_end) {
  seqtab <- makeSequenceTable(dada_F)
} else {
  dada_R <- dada(derep_R, err=errR, multithread=TRUE)
  
  # Merge paired reads
  mergers <- mergePairs(dada_F, derep_F, dada_R, derep_R, verbose=TRUE)
  seqtab <- makeSequenceTable(mergers)
}

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save ASV table
saveRDS(seqtab.nochim, file=output_asv)

# ----------------
# RDP Taxonomic Classification
# ----------------

# Provide path to reference database
rdp_trainset <- "rdp_train_set_19.fa.gz"

if (file.exists(rdp_trainset)) {
  tax_table <- assignTaxonomy(seqtab.nochim, rdp_trainset, multithread=TRUE)
  saveRDS(tax_table, file=output_tax)
  cat("RDP classification completed. Results saved in", output_tax, "\n")
} else {
  cat("RDP database not found! Skipping taxonomy assignment.\n")
}

cat("DADA2 processing completed. ASV table saved in", output_asv, "\n")
