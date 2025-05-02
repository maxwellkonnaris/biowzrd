#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(BiocParallel)
})

set.seed(100)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: run_dada2.R <bioproject> <inputfile> <threads> <outdir>")

biop       <- args[1]
inputfile  <- args[2]
outdir     <- args[4]
chkdir     <- file.path(outdir, "checkpoints")
dir.create(chkdir, showWarnings=FALSE, recursive=TRUE)
threads <- as.integer(args[3])
register(MulticoreParam(workers = threads))
options(mc.cores = threads) 
VERBOSE <- TRUE

manifest <- read.csv(inputfile, stringsAsFactors = FALSE)
sub      <- subset(manifest, Bioproject == biop)
if (nrow(sub) == 0) {
  stop("No rows found in input file for Bioproject: ", biop)
}

fnF_list <- sub$Fastq1
if (any(!file.exists(fnF_list))) {
  missing <- fnF_list[!file.exists(fnF_list)]
  stop("The following FASTQ files are missing:\n", paste(missing, collapse = "\n"))
}

message("------------------------- Start Bioproject:", biop)

# 1) FILTER & TRIM -------------------------------------------------------
filtpath <- file.path(outdir, "filtered")
dir.create(filtpath, showWarnings = FALSE, recursive = TRUE)

SIZE_THRESHOLD_GB <- 20  # Chunk if larger than this
CHUNK_LINES <- 4e6       # Number of lines per chunk (1M reads = 4M lines)

for (fnF in fnF_list) {
  if (!file.exists(fnF)) {
    message("File does not exist: ", fnF)
    next
  }

  base <- sub("\\.fastq\\.gz$", "", basename(fnF))
  outF <- file.path(filtpath, paste0("filtered_", base, ".fastq"))
  chk  <- file.path(chkdir, paste0("filter_", base, ".rds"))

  if (file.exists(chk)) {
    message("Skipping filter: ", basename(fnF), " (already done)")
    next
  }

  fsize_gb <- file.info(fnF)$size / 1e9
  message(sprintf("Filtering: %s (size=%.1f GB)", basename(fnF), fsize_gb))

  # Large file: chunk it first
  if (fsize_gb > SIZE_THRESHOLD_GB) {
    tmpdir <- tempfile(pattern = paste0("chunks_", base))
    dir.create(tmpdir)
    message("Chunking large file with seqkit...")

    # Run seqkit to split file into chunks of ~1M reads
    cmd <- sprintf("seqkit split -s 1e6 -O '%s' '%s'", tmpdir, fnF)
    system(cmd)

    chunk_files <- list.files(tmpdir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
    filtered_chunks <- c()

    for (chunk in chunk_files) {
      out_chunk <- tempfile(tmpdir = tmpdir, fileext = ".fastq")
      message("  ↳ Filtering chunk: ", basename(chunk))
      filterAndTrim(
        chunk, out_chunk,
        truncLen    = 0,
        maxEE       = 10,
        truncQ      = 2,
        maxN        = 0,
        rm.phix     = TRUE,
        compress    = FALSE,
        multithread = TRUE,
        verbose     = VERBOSE
      )
      filtered_chunks <- c(filtered_chunks, out_chunk)
    }

    # Merge filtered chunks
    message("Merging filtered chunks into final output: ", outF)
    file.create(outF)
    for (fc in filtered_chunks) {
      writeLines(readLines(fc), outF, sep = "\n", useBytes = TRUE, append = TRUE)
    }

    unlink(tmpdir, recursive = TRUE)

  } else {
    # Small enough, filter directly
    filterAndTrim(
      fnF, outF,
      truncLen    = 0,
      maxEE       = 10,
      truncQ      = 2,
      maxN        = 0,
      rm.phix     = TRUE,
      compress    = FALSE,
      multithread = TRUE,
      verbose     = VERBOSE
    )
  }

  saveRDS(TRUE, chk)
}

message("------------------------- Finished Filtering:", biop)

# 2) LEARN ERRORS --------------------------------------------------------

err_chk <- file.path(chkdir, "errors.rds")
smp_chk <- file.path(chkdir, "sample_names.rds")

# Find filtered FASTQ files (already generated or new ones)
filts <- list.files(filtpath, pattern="filtered_.*\\.fastq$", full.names=TRUE)
sample.names <- sub(
  "^filtered_([^_]+)(?:_1)?\\.fastq$",
  "\\1",
  basename(filts),
  perl = TRUE
)

# Check if an existing error model is out of date
if (file.exists(smp_chk)) {
  saved_sample.names <- readRDS(smp_chk)
  
  if (!identical(sort(saved_sample.names), sort(sample.names))) {
    message("Sample list has changed since last run. Clearing old error model.")
    unlink(c(err_chk, smp_chk))
  }
}

# Reload or re-learn the error model
if (file.exists(err_chk) && file.exists(smp_chk)) {
  message("Loading existing error model and sample names")
  err <- readRDS(err_chk)
  sample.names <- readRDS(smp_chk)
  filts <- list.files(filtpath, pattern="filtered_.*\\.fastq$", full.names=TRUE)
  names(filts) <- sample.names
  if (anyDuplicated(sample.names)) {
    stop("Duplicate sample names detected!")
  }
} else {
  message("Learning error rates")
  filts <- list.files(filtpath, pattern="filtered_.*\\.fastq$", full.names=TRUE)
  sample.names <- sub(
    "^filtered_([^_]+)(?:_1)?\\.fastq$",
    "\\1",
    basename(filts),
    perl = TRUE
  )
  names(filts) <- sample.names
  if (anyDuplicated(sample.names)) {
    stop("Duplicate sample names detected!")
  }

  err <- learnErrors(filts, nbases=1e8, multithread=TRUE, randomize=TRUE, verbose=VERBOSE)
  saveRDS(err, err_chk)
  saveRDS(sample.names, smp_chk)
}
message("------------------------- Finished Error Model: ", biop)

# 3) dada() PER-SAMPLE ----------------------------------------------------

dds_chk <- file.path(chkdir, "dds_list.rds")
if (file.exists(dds_chk)) {
  message("Loading existing DADA2 results")
  dds <- readRDS(dds_chk)
} else {
  message("Running DADA2 per sample")

  dds <- vector("list", length(sample.names))
  names(dds) <- sample.names

  for (sam in sample.names) {
    cat("→ Processing:", sam, "\n")
    derep      <- derepFastq(filts[[sam]], verbose=VERBOSE)
    dds[[sam]] <- dada(
      derep,
      err         = err,
      multithread = TRUE,
      verbose     = VERBOSE
    )
  }

  saveRDS(dds, dds_chk)
}

message("------------------------- Finished Derep and Dada:", biop)

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

message("------------------------- Finished Sequence Table:", biop)

# 5) REMOVE CHIMERAS ------------------------------------------------------

nochim_chk <- file.path(chkdir, "seqtab_nochim.rds")
if (file.exists(nochim_chk)) {
  message("Loading chimera-free table")
  seqtab_nochim <- readRDS(nochim_chk)
} else {
  message("Removing chimeras")
  seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose = VERBOSE)
  saveRDS(seqtab_nochim, nochim_chk)
}

message("------------------------- Finished Chimeras:", biop)

# 6) TAXONOMY -------------------------------------------------------------

if (!file.exists("rdp_19_toGenus_trainset.fa.gz")) {
    stop("RDP taxonomy file not found. Please download or specify the correct path.")
}

tax_chk <- file.path(chkdir, "tax.rds")
if (file.exists(tax_chk)) {
  message("Loading taxonomy")
  tax <- readRDS(tax_chk)
} else {
  message("Assigning taxonomy")
  tax <- assignTaxonomy(seqtab_nochim, "rdp_19_toGenus_trainset.fa.gz", multithread=TRUE, verbose = VERBOSE)
  saveRDS(tax, tax_chk)
}

message("------------------------- Finished Taxonomy:", biop)

# 7) FINAL SAVES ----------------------------------------------------------

saveRDS(seqtab,        file.path(outdir, paste0("asv_", biop, ".rds")))
saveRDS(seqtab_nochim, file.path(outdir, paste0(biop, "_dada2_counts.rds")))
saveRDS(tax,           file.path(outdir, paste0(biop, "_dada2_taxa.rds")))

message("All steps complete.")
quit(save="no", status=0)
