#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript merge_dada2.R <bioproject> <seqtab_file1> [<seqtab_file2> ...]")
}

bioproject   <- args[1]
seqtab_files <- args[-1]
rdp_trainset <- "rdp_19_toGenus_trainset.fa.gz"

# ---- validation -----------------------------------------------------------
if (!file.exists(rdp_trainset)) stop("RDP database not found: ", rdp_trainset)
if (any(!file.exists(seqtab_files))) stop("One or more sequence table files do not exist.")

# ---- merge + chimera removal ----------------------------------------------
seqtabs       <- lapply(seqtab_files, readRDS)
seqtab_merged <- do.call(mergeSequenceTables, seqtabs)
seqtab_nochim <- removeBimeraDenovo(seqtab_merged,
                                    method = "consensus",
                                    minFoldParentOverAbundance = 3.5,
                                    multithread = 4, verbose = TRUE)

# ---- save merged ASV table -----------------------------------------------
out_dir           <- dirname(seqtab_files[1])
merged_asv_file   <- file.path(out_dir, paste0(bioproject, "_dada2_merged_nochim.rds"))
saveRDS(seqtab_nochim, merged_asv_file)
cat("Saved merged ASV table to:", merged_asv_file, "\n")

# ---- one-time global taxonomy assignment ---------------------------------
asv_seqs           <- colnames(seqtab_nochim)
tax_table_merged   <- assignTaxonomy(asv_seqs,
                                     rdp_trainset,
                                     multithread = 4,
                                     verbose = TRUE)

merged_tax_file    <- file.path(out_dir, paste0(bioproject, "_dada2_taxonomy_merged.rds"))
saveRDS(tax_table_merged, merged_tax_file)
cat("Saved merged taxonomy table to:", merged_tax_file, "\n")

cat("Merge, chimera removal, and global taxonomy assignment completed for", bioproject, "\n")
