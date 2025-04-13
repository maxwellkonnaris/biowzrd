#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_dada2.R <forward_fastq> [<reverse_fastq>|--seqtab <seqtab.rds>] [--output-asv <asv.rds>] [--output-tax <tax.rds>]", call.=FALSE)
}

# Parse arguments
mode <- "process"  # Default: process FASTQs
input_F <- args[1]
input_R <- NULL
seqtab_file <- NULL
output_asv <- NULL
output_tax <- NULL

if (args[2] == "--seqtab") {
  mode <- "taxonomy"
  seqtab_file <- args[3]
  for (i in 4:length(args)) {
    if (args[i] == "--output-asv") output_asv <- args[i+1]
    if (args[i] == "--output-tax") output_tax <- args[i+1]
  }
} else {
  input_R <- args[2]
}

# Single-end or paired-end
single_end <- is.null(input_R) && mode == "process"

if (mode == "process") {
  # Output filenames
  sample_name <- gsub("_1\\.fastq\\.gz$|.fastq.gz$", "", basename(input_F))
  output_seqtab <- paste0("seqtab_", sample_name, ".rds")

  # Define paths
  filtered_path <- "filtered_fastq/"
  if (!dir.exists(filtered_path)) dir.create(filtered_path)

  filt_F <- file.path(filtered_path, paste0("filtered_", basename(input_F)))
  filt_R <- if (!single_end) file.path(filtered_path, paste0("filtered_", basename(input_R))) else NULL

  # Filtering and trimming
  if (single_end) {
    out <- filterAndTrim(input_F, filt_F, truncLen=250, maxN=0, maxEE=2, truncQ=2, 
                         rm.phix=TRUE, compress=TRUE, multithread=TRUE)
  } else {
    out <- filterAndTrim(fwd=input_F, filt=filt_F, rev=input_R, filt.rev=filt_R, 
                         truncLen=c(250,200), maxN=0, maxEE=c(2,2), truncQ=2, 
                         rm.phix=TRUE, compress=TRUE, multithread=TRUE)
  }

  # Check if any reads passed filtering
  if (sum(out[, "reads.out"]) == 0) {
    stop("No reads passed quality filtering for ", sample_name)
  }

  # Learn error rates
  errF <- learnErrors(filt_F, multithread=TRUE)
  if (!single_end) {
    errR <- learnErrors(filt_R, multithread=TRUE)
  }

  # Dereplication
  derep_F <- derepFastq(filt_F, verbose=TRUE)
  names(derep_F) <- sample_name

  if (!single_end) {
    derep_R <- derepFastq(filt_R, verbose=TRUE)
    names(derep_R) <- sample_name
  }

  # DADA2 core processing
  dada_F <- dada(derep_F, err=errF, multithread=TRUE)

  if (single_end) {
    seqtab <- makeSequenceTable(dada_F)
  } else {
    dada_R <- dada(derep_R, err=errR, multithread=TRUE)
    mergers <- mergePairs(dada_F, derep_F, dada_R, derep_R, verbose=TRUE)
    seqtab <- makeSequenceTable(mergers)
  }

  # Save sequence table (pre-chimera removal)
  saveRDS(seqtab, file=output_seqtab)
  cat("Sequence table saved in", output_seqtab, "\n")

  # Clean up filtered FASTQs
  unlink(filtered_path, recursive=TRUE)
} else if (mode == "taxonomy") {
  # Taxonomy assignment
  if (is.null(output_asv) || is.null(output_tax)) {
    stop("Both --output-asv and --output-tax must be specified in taxonomy mode")
  }

  seqtab <- readRDS(seqtab_file)
  rdp_trainset <- "rdp_19_toGenus_trainset.fa.gz"

  if (file.exists(rdp_trainset)) {
    tax_table <- assignTaxonomy(seqtab, rdp_trainset, multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab, file=output_asv)
    saveRDS(tax_table, file=output_tax)
    cat("Taxonomy assignment completed. ASV table saved in", output_asv, ", taxonomy in", output_tax, "\n")
  } else {
    cat("RDP database not found! Skipping taxonomy assignment.\n")
  }
}
