#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(BiocParallel)
})

set.seed(100)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: run_dada2_partial.R <bioproject> <comma,separated,forward.fastq.gz> <threads>")

biop       <- args[1]
fnF_list   <- strsplit(args[2], ",")[[1]]
chkdir     <- file.path(biop, "checkpoints")
dir.create(chkdir, showWarnings=FALSE, recursive=TRUE)
threads <- as.integer(Sys.getenv("D2_THREADS", 4))
register(MulticoreParam(workers = threads))

# 1) FILTER & TRIM -------------------------------------------------------

filtpath <- file.path(biop, "filtered")
dir.create(filtpath, showWarnings=FALSE, recursive=TRUE)

for(fnF in fnF_list) {
  outF <- file.path(filtpath, paste0("filtered_", basename(fnF)))
  chk  <- file.path(chkdir, paste0("filter_", basename(fnF), ".rds"))
  
  if (file.exists(chk)) {
    message("Skipping filter: ", basename(fnF), " (already done)")
  } else {
    message("Filtering: ", basename(fnF))
    filterAndTrim(fnF, outF,
                  truncLen = 0, maxEE = 10, truncQ = 2, maxN = 0,
                  rm.phix = TRUE, compress = FALSE,
                  multithread = TRUE, verbose = TRUE)
    saveRDS(TRUE, chk)  # simply mark this file as done
  }
}

# 2) LEARN ERRORS --------------------------------------------------------

err_chk <- file.path(chkdir, "errors.rds")
if (file.exists(err_chk)) {
  message("Loading existing error model")
  err <- readRDS(err_chk)
} else {
  message("Learning error rates")
  filts       <- list.files(filtpath, pattern="fastq", full.names=TRUE)
  sample.names<- sapply(strsplit(basename(filts), "_"), `[`, 1)
  names(filts) <- sample.names
  err <- learnErrors(filts, nbases=1e8, multithread=threads, randomize=TRUE)
  saveRDS(err, err_chk)
}

# 3) dada() PER-SAMPLE ----------------------------------------------------

dds_chk <- file.path(chkdir, "dds_list.rds")
if (file.exists(dds_chk)) {
  message("Loading existing DADA2 results")
  dds <- readRDS(dds_chk)
} else {
  message("Running DADA2 per sample")
  dds <- vector("list", length(names(err)))
  names(dds) <- names(err)
  for (sam in names(err)) {
    message("  → ", sam)
    derep      <- derepFastq(filts[[sam]])
    dds[[sam]]<- dada(derep, err=err, multithread=threads)
  }
  saveRDS(dds, dds_chk)
}

# 4) MAKE SEQUENCE TABLE --------------------------------------------------

seqtab_chk <- file.path(chkdir, "seqtab.rds")
if (file.exists(seqtab_chk)) {
  message("Loading existing sequence table")
  seqtab <- readRDS(seqtab_chk)
} else {
  message("Building sequence table")
  seqtab <- makeSequenceTable(dds)
  saveRDS(seqtab, seqtab_chk)
}

# 5) REMOVE CHIMERAS ------------------------------------------------------

nochim_chk <- file.path(chkdir, "seqtab_nochim.rds")
if (file.exists(nochim_chk)) {
  message("Loading chimera-free table")
  seqtab_nochim <- readRDS(nochim_chk)
} else {
  message("Removing chimeras")
  seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads)
  saveRDS(seqtab_nochim, nochim_chk)
}

# 6) TAXONOMY -------------------------------------------------------------

tax_chk <- file.path(chkdir, "tax.rds")
if (file.exists(tax_chk)) {
  message("Loading taxonomy")
  tax <- readRDS(tax_chk)
} else {
  message("Assigning taxonomy")
  tax <- assignTaxonomy(seqtab_nochim, "rdp_19_toGenus_trainset.fa.gz", multithread=threads)
  saveRDS(tax, tax_chk)
}

# 7) FINAL SAVES ----------------------------------------------------------

saveRDS(seqtab,        file.path(biop, paste0("asv_", biop, ".rds")))
saveRDS(seqtab_nochim, file.path(biop, paste0(biop, "_dada2_counts.rds")))
saveRDS(tax,           file.path(biop, paste0(biop, "_dada2_taxa.rds")))

message("All steps complete.")
quit(save="no", status=0)
