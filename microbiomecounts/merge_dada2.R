#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript merge_dada2.R <bioproject> <seqtab_file1> [<seqtab_file2> ...]", call.=FALSE)
}

bioproject <- args[1]
seqtab_files <- args[-1]

# Load and merge sequence tables
seqtabs <- lapply(seqtab_files, readRDS)
seqtab_merged <- do.call(mergeSequenceTables, seqtabs)

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab_merged, method="consensus", minFoldParentOverAbundance=3.5, multithread=TRUE, verbose=TRUE)

# Split back to per-accession tables and assign taxonomy
sample_names <- sub("^seqtab_([^.]+)\\.rds$", "\\1", basename(seqtab_files))
rdp_trainset <- "rdp_19_toGenus_trainset.fa.gz"

for (i in seq_along(sample_names)) {
  sample <- sample_names[i]
  seqtab_sample <- seqtab_nochim[, , sample, drop=FALSE]
  output_asv <- file.path(bioproject, paste0("dada2_results_", sample, ".rds"))
  output_tax <- file.path(bioproject, paste0("dada2_taxonomy_", sample, ".rds"))
  
  if (file.exists(rdp_trainset)) {
    tax_table <- assignTaxonomy(seqtab_sample, rdp_trainset, multithread=TRUE, verbose=TRUE)
    saveRDS(seqtab_sample, file=output_asv)
    saveRDS(tax_table, file=output_tax)
    cat("Processed", sample, ": ASV table saved in", output_asv, ", taxonomy in", output_tax, "\n")
  } else {
    cat("RDP database not found! Skipping taxonomy for", sample, "\n")
    saveRDS(seqtab_sample, file=output_asv)
  }
}

cat("Chimera removal and taxonomy assignment completed for", bioproject, "\n")
