#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: Rscript run_dada2_partial.R <forward_fastq> [<reverse_fastq>] <output_file> [<filter_mode>]")
}

parse_args <- function(a) {
  n <- length(a)
  if (n == 2) {
    list(fnF = a[1], fnR = "", out = a[2], mode = "light")
  } else if (n == 3 && grepl("\\.fastq", a[2], ignore.case = TRUE)) {
    list(fnF = a[1], fnR = a[2], out = a[3], mode = "light")
  } else if (n == 3) {
    list(fnF = a[1], fnR = "", out = a[2], mode = tolower(a[3]))
  } else if (n == 4) {
    list(fnF = a[1], fnR = a[2], out = a[3], mode = tolower(a[4]))
  } else stop("Unrecognised argument pattern")
}

p <- parse_args(args)
fnF <- p$fnF
fnR <- p$fnR
output_file <- p$out
filter_mode <- p$mode

if (!filter_mode %in% c("strict", "light", "none")) filter_mode <- "light"

sample <- sub("^asv_([^.]+)\\.rds$", "\\1", basename(output_file))
out_dir <- dirname(output_file)

if (!file.exists(fnF)) stop("Forward FASTQ does not exist: ", fnF)
if (nzchar(fnR) && !file.exists(fnR)) stop("Reverse FASTQ does not exist: ", fnR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

filt_F <- file.path(out_dir, paste0("filtered_", basename(fnF)))
filt_R <- if (nzchar(fnR)) file.path(out_dir, paste0("filtered_", basename(fnR))) else NULL
single_end <- TRUE

if (filter_mode == "strict") {
  truncLen <- 200; maxEE <- 2; truncQ <- 2
} else if (filter_mode == "light") {
  truncLen <- 0; maxEE <- 10; truncQ <- 2
} else {
  truncLen <- 0; maxEE <- Inf; truncQ <- 0
}

timings <- list()

t0 <- proc.time()
out <- filterAndTrim(fnF, filt_F,
                     truncLen = truncLen, maxN = 0, maxEE = maxEE, truncQ = truncQ,
                     rm.phix = TRUE, compress = FALSE, multithread = 4)
timings$filterAndTrim <- (proc.time() - t0)[["elapsed"]]

if (sum(out[, "reads.out"]) == 0) stop("No reads passed quality filtering for ", sample)

t0 <- proc.time()
derepForErr <- derepFastq(filt_F, verbose = FALSE)
if (sum(derepForErr$quals[[1]] != "") > 1e6) {
  set.seed(42)
  sample_ix <- sample(seq_along(derepForErr$quals), size = 1e6)
  derepForErr <- list(quals = derepForErr$quals[sample_ix])
}
errF <- learnErrors(derepForErr, multithread = 4)
timings$learnErrors <- (proc.time() - t0)[["elapsed"]]

t0 <- proc.time()
derep_F <- derepFastq(filt_F)
timings$derepFastq <- (proc.time() - t0)[["elapsed"]]

t0 <- proc.time()
dada_F <- suppressMessages(dada(derep_F, err = errF, multithread = 4))
timings$dada <- (proc.time() - t0)[["elapsed"]]

t0 <- proc.time()
seqtab <- makeSequenceTable(dada_F)
timings$makeSequenceTable <- (proc.time() - t0)[["elapsed"]]

if (ncol(seqtab) > 0) {
  message(sprintf("[DADA2] %s → Single-end: %d ASVs", sample, ncol(seqtab)))
  
  t0 <- proc.time()
  saveRDS(seqtab, output_file)
  write.table(out, file.path(out_dir, paste0(sample, "_filter_summary.txt")),
              sep = "\t", quote = FALSE)
  unlink(filt_F)
  timings$saveOutputs <- (proc.time() - t0)[["elapsed"]]
  
  message(sprintf(
    "[%s] TIMING SUMMARY: filterAndTrim=%.1fs, learnErrors=%.1fs, derepFastq=%.1fs, dada=%.1fs, makeSequenceTable=%.1fs, saveOutputs=%.1fs",
    sample,
    timings$filterAndTrim, timings$learnErrors, timings$derepFastq,
    timings$dada, timings$makeSequenceTable, timings$saveOutputs
  ))
  
  quit(save = "no", status = 0)
} else {
  stop(sprintf("[DADA2] %s → No ASVs detected (even in fallback single-end)", sample))
}


# ---- OPTIONAL: Paired-end fallback (currently commented out) ----
# if (!single_end) {
#   out <- filterAndTrim(fnF, filt_F, fnR, filt_R,
#                        truncLen = c(truncLen, truncLen), maxN = 0,
#                        maxEE = maxEE, truncQ = truncQ,
#                        rm.phix = TRUE, compress = FALSE, multithread = 4)
#   if (sum(out[, "reads.out"]) == 0) stop("No reads passed quality filtering for ", sample)
#   errR <- learnErrors(filt_R, nbases=1e8, multithread = 4)
#   derep_R <- derepFastq(filt_R)
#   dada_R <- suppressMessages(dada(derep_R, err = errR, multithread = 4))
#   mergers <- suppressMessages(mergePairs(dada_F, derep_F, dada_R, derep_R, verbose = FALSE))
#   seqtab <- makeSequenceTable(mergers)
#   if (ncol(seqtab) > 0) {
#     message(sprintf("[DADA2] %s → Paired-end: %d merged ASVs", sample, ncol(seqtab)))
#     saveRDS(seqtab, output_file)
#     write.table(out, file.path(out_dir, paste0(sample, "_filter_summary.txt")),
#                 sep = "\t", quote = FALSE)
#     unlink(c(filt_F, filt_R))
#     quit(save = "no", status = 0)
#   } else {
#     message(sprintf("[DADA2] %s → Paired-end failed, falling back to single-end", sample))
#   }
# }
