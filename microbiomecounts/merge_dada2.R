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
if (!file.exists(rdp_trainset))           stop("RDP database not found: ", rdp_trainset)
if (any(!file.exists(seqtab_files)))      stop("One or more sequence table files do not exist.")

# ---- merge + chimera removal ----------------------------------------------
seqtabs        <- lapply(seqtab_files, \(f) readRDS(f))
seqtab_merged  <- do.call(mergeSequenceTables, seqtabs)
seqtab_nochim  <- removeBimeraDenovo(seqtab_merged, method = "consensus",
                                     minFoldParentOverAbundance = 3.5,
                                     multithread = TRUE, verbose = TRUE)

out_dir        <- dirname(seqtab_files[1])
saveRDS(seqtab_nochim,
        file.path(out_dir, paste0(bioproject, "_dada2_merged_nochim.rds")))

# ---- per‑sample taxonomy --------------------------------------------------
sample_names <- sub("^asv_([^.]+)\\.rds$", "\\1", basename(seqtab_files))

for (i in seq_along(sample_names)) {
  sample_id      <- sample_names[i]
  seqtab_sample  <- seqtab_nochim[sample_id, , drop = FALSE]

  if (nrow(seqtab_sample) == 0) {
    warning("No reads for sample ", sample_id, "; skipping taxonomy.")
    next
  }

  asv_file <- file.path(out_dir, paste0("dada2_results_",   sample_id, ".rds"))
  tax_file <- file.path(out_dir, paste0("dada2_taxonomy_", sample_id, ".rds"))

  tax_table <- tryCatch(
    assignTaxonomy(seqtab_sample, rdp_trainset, multithread = TRUE, verbose = TRUE),
    error = \(e) { warning("Taxonomy failed for ", sample_id, ": ", conditionMessage(e)); NULL }
  )

  saveRDS(seqtab_sample, asv_file)
  if (!is.null(tax_table)) saveRDS(tax_table, tax_file)

  cat("Processed", sample_id,
      "→ ASV:", basename(asv_file),
      if (!is.null(tax_table)) paste0(" taxonomy:", basename(tax_file)) else "(no taxonomy)",
      "\n")
}

cat("Merge, chimera removal, and taxonomy assignment completed for", bioproject, "\n")
