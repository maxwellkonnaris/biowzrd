#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_dada2.R <input_fastq>", call.=FALSE)
}

input_file <- args[1]
output_file <- paste0("dada2_results_", gsub(".fastq.gz", "", basename(input_file)), ".rds")

# Define paths
filt_file <- paste0("filtered_", basename(input_file))
filtered_path <- "filtered_fastq/"

# Create directory if not exists
if (!dir.exists(filtered_path)) dir.create(filtered_path)

# Filtering and trimming
filt_path <- file.path(filtered_path, filt_file)
out <- filterAndTrim(input_file, filt_path, truncLen=c(250,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Learn error rates
errF <- learnErrors(filt_path, multithread=TRUE)
errR <- learnErrors(filt_path, multithread=TRUE)

# Dereplicate
derep <- derepFastq(filt_path)
names(derep) <- gsub(".fastq.gz", "", basename(filt_path))

# DADA2 core processing
dadaFs <- dada(derep, err=errF, multithread=TRUE)
seqtab <- makeSequenceTable(dadaFs)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save results
saveRDS(seqtab.nochim, file=output_file)

cat("DADA2 processing completed. Results saved in", output_file, "\n")
