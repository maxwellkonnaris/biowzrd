#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: Rscript run_dada2_partial.R <forward_fastq> [<reverse_fastq>] <output_file> [<filter_mode>]")
}

# --------------------------------------------------------------------------
# Robust argument parsing
# --------------------------------------------------------------------------
parse_args <- function(a) {
  n <- length(a)
  if (n == 2) {                             # single‑end, default mode
    list(fnF = a[1], fnR = "",  out = a[2], mode = "light")
  } else if (n == 3 && grepl("\\.fastq", a[2], ignore.case = TRUE)) {
    list(fnF = a[1], fnR = a[2], out = a[3], mode = "light")   # paired‑end
  } else if (n == 3) {                      # single‑end + mode
    list(fnF = a[1], fnR = "",  out = a[2], mode = tolower(a[3]))
  } else if (n == 4) {                      # paired‑end + mode
    list(fnF = a[1], fnR = a[2], out = a[3], mode = tolower(a[4]))
  } else stop("Unrecognised argument pattern")
}

p          <- parse_args(args)
fnF        <- p$fnF
fnR        <- p$fnR
output_file<- p$out
filter_mode<- p$mode

if (!filter_mode %in% c("strict", "light", "none")) filter_mode <- "light"

sample   <- sub("^asv_([^.]+)\\.rds$", "\\1", basename(output_file))
out_dir  <- dirname(output_file)

# ---- I/O checks -----------------------------------------------------------
if (!file.exists(fnF)) stop("Forward FASTQ does not exist: ", fnF)
if (nzchar(fnR) && !file.exists(fnR)) stop("Reverse FASTQ does not exist: ", fnR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

filt_F <- file.path(out_dir, paste0("filtered_", basename(fnF)))
filt_R <- if (nzchar(fnR)) file.path(out_dir, paste0("filtered_", basename(fnR))) else NULL
single_end <- !nzchar(fnR)

# ---- mode‑specific parameters --------------------------------------------
if (filter_mode == "strict") {
  truncLen <- if (single_end) 200 else c(200, 160); maxEE <- 2;   truncQ <- 2
} else if (filter_mode == "light") {
  truncLen <- if (single_end)   0 else c(0, 0);   maxEE <- 10;   truncQ <- 2
} else {                                   # "none"
  truncLen <- if (single_end)   0 else c(0, 0);   maxEE <- Inf;  truncQ <- 0
}

# ---- filter & trim --------------------------------------------------------
message("Filtering (mode = ", filter_mode, ")")
out <- if (single_end) {
  filterAndTrim(fnF, filt_F,               truncLen = truncLen,
                maxN = 0, maxEE = maxEE, truncQ = truncQ,
                rm.phix = TRUE, compress = FALSE, multithread = TRUE)
} else {
  filterAndTrim(fnF, filt_F, fnR, filt_R, truncLen = truncLen,
                maxN = 0, maxEE = maxEE, truncQ = truncQ,
                rm.phix = TRUE, compress = FALSE, multithread = TRUE)
}

if (sum(out[, "reads.out"]) == 0) stop("No reads passed quality filtering for ", sample)

# ---- learn errors --------------------------------------------------------
errF <- learnErrors(filt_F, multithread = TRUE)
if (!single_end) errR <- learnErrors(filt_R, multithread = TRUE)

# ---- dereplication --------------------------------------------------------
derep_F <- derepFastq(filt_F); names(derep_F) <- sample
if (!single_end) { derep_R <- derepFastq(filt_R); names(derep_R) <- sample }

# ---- DADA inference & table formation ------------------------------------
dada_F <- dada(derep_F, err = errF, multithread = TRUE)

seqtab <- if (single_end) {
  makeSequenceTable(dada_F)
} else {
  dada_R  <- dada(derep_R, err = errR, multithread = TRUE)
  mergers <- mergePairs(dada_F, derep_F, dada_R, derep_R, verbose = TRUE)
  makeSequenceTable(mergers)
}

# ---- save & clean ---------------------------------------------------------
saveRDS(seqtab, output_file)
write.table(out, file.path(out_dir, paste0(sample, "_filter_summary.txt")),
            sep = "\t", quote = FALSE)

unlink(filt_F); if (!single_end) unlink(filt_R)
