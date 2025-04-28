#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dada2)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
usage <- "
Usage:
  Rscript run_dada2_partial.R <forward.fastq.gz> [<reverse.fastq.gz>] <out.rds> [strict|light|none]
"
if (length(args) < 2 || length(args) > 4) stop(usage)

# Parse arguments
if (grepl("\\.fastq", args[2], ignore.case = TRUE)) {
  fnF <- args[1]; fnR <- args[2]; out <- args[3]
  mode <- if (length(args) == 4) tolower(args[4]) else "light"
} else {
  fnF <- args[1]; fnR <- ""
  out <- args[2]; mode <- if (length(args) == 3) tolower(args[3]) else "light"
}
stopifnot(mode %in% c("strict", "light", "none"))

# Filtering parameters
pars <- switch(
  mode,
  strict = list(truncLen = 200, maxEE = 2,  truncQ = 2),
  light  = list(truncLen = 0,   maxEE = 10, truncQ = 2),
  list(truncLen = 0, maxEE = Inf, truncQ = 0)
)

fnR <- ""

# Setup parallelization
threads <- as.integer(Sys.getenv("D2_THREADS", 4))
register(MulticoreParam(workers = threads))
message(sprintf("[DADA2] Using %d threads (via BiocParallel)", threads))

# Prepare paths and sample name
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
filt_F <- file.path(dirname(out), paste0("filtered_", basename(fnF)))
filt_R <- if (nzchar(fnR)) file.path(dirname(out), paste0("filtered_", basename(fnR))) else NULL
sample <- sub("^asv_([^.]+)\\.rds$", "\\1", basename(out))

timings <- list()

# 1. filterAndTrim
t0 <- proc.time()
out_filt <- filterAndTrim(
  fnF, filt_F,
  if (nzchar(fnR)) fnR else NULL,
  if (nzchar(fnR)) filt_R else NULL,
  truncLen = pars$truncLen,
  maxN = 0,
  maxEE = pars$maxEE,
  truncQ = pars$truncQ,
  rm.phix = TRUE,
  compress = FALSE,
  multithread = FALSE
)
timings$filterAndTrim <- (proc.time() - t0)[["elapsed"]]
if (sum(out_filt[, "reads.out"]) == 0) stop("No reads passed filtering for ", sample)

# 3. derepFastq
t0 <- proc.time()
derep_F <- derepFastq(filt_F)
if (nzchar(fnR)) derep_R <- derepFastq(filt_R)
timings$derepFastq <- (proc.time() - t0)[["elapsed"]]

# 2. learnErrors
# Single-end or paired: learn on each filtered file
learn_errors <- function(fq) {
  t1 <- proc.time()
  err <- learnErrors(fq, nbases = 1e8, randomize = TRUE, multithread=threads)
  elapsed <- (proc.time() - t1)[["elapsed"]]
  list(err = err, time = elapsed)
}

leF <- learn_errors(derep_F)
timings$learnErrors_F <- leF$time

if (nzchar(fnR)) {
  leR <- learn_errors(derep_R)
  timings$learnErrors_R <- leR$time
}

# 4. dada
t0 <- proc.time()
dada_F <- dada(derep_F, err = leF$err, multithread=threads)
if (nzchar(fnR)) dada_R <- dada(derep_R, err = leR$err, multithread=threads)
timings$dada <- (proc.time() - t0)[["elapsed"]]

# 5. mergePairs or single-end
t0 <- proc.time()
if (nzchar(fnR)) {
  mergers <- mergePairs(dada_F, derep_F, dada_R, derep_R, verbose = FALSE)
  seqtab <- makeSequenceTable(mergers)
  timings$mergePairs <- (proc.time() - t0)[["elapsed"]]
} else {
  seqtab <- makeSequenceTable(dada_F)
  timings$makeSequenceTable <- (proc.time() - t0)[["elapsed"]]
}

# 6. Save and cleanup
if (ncol(seqtab) == 0) stop(sprintf("[DADA2] %s → No ASVs detected", sample))
t0 <- proc.time()
saveRDS(seqtab, out)
write.table(out_filt,
            file.path(dirname(out), paste0(sample, "_filter_summary.txt")),
            sep = "\t", quote = FALSE)
unlink(c(filt_F, filt_R))
timings$saveOutputs <- (proc.time() - t0)[["elapsed"]]

# 7. Report timings
msg <- sprintf(
  "[%s] TIMING: filter=%.1fs, errF=%.1fs%s, derep=%.1fs, dada=%.1fs%s%s save=%.1fs", 
  sample,
  timings$filterAndTrim,
  timings$learnErrors_F,
  if (nzchar(fnR)) sprintf(", errR=%.1fs", timings$learnErrors_R) else "",
  timings$derepFastq,
  timings$dada,
  if (nzchar(fnR)) sprintf(", merge=%.1fs, seqtab=%0.1fs", timings$mergePairs, timings$makeSequenceTable) else "",
  timings$saveOutputs
)
message(msg)

quit(save = "no", status = 0)



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
