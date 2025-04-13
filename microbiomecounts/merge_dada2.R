#!/usr/bin/env Rscript

library(dada2)

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
seqtab_nochim <- removeBimeraDenovo(seqtab_merged, method="consensus", multithread=TRUE, verbose=TRUE)

# Split back to per-accession tables
sample_names <- sub("^seqtab_([^.]+)\\.rds$", "\\1", basename(seqtab_files))
for (i in seq_along(sample_names)) {
  sample <- sample_names[i]
  seqtab_sample <- seqtab_nochim[, , sample, drop=FALSE]
  saveRDS(seqtab_sample, file.path(bioproject, paste0("seqtab_nochim_", sample, ".rds")))
}

cat("Chimera removal completed for", bioproject, "\n")
