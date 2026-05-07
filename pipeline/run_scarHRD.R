#!/usr/bin/env Rscript
# run_scarHRD.R — compute HRD scores via scarHRD for all input files in a directory.
#
# Usage:
#   Rscript run_scarHRD.R <input_dir> <output_dir>
#
# Reads all *scarHRD_input*.txt files from <input_dir>, runs scar_score() on each,
# and writes ALL_scarHRD_scores.csv to <output_dir>.
#
# scarHRD creates intermediate files in the working directory, so each sample
# is processed in a temporary subdirectory that is cleaned up afterwards.

suppressPackageStartupMessages(library(scarHRD))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript run_scarHRD.R <input_dir> <output_dir>\n")
  quit(status = 1)
}

input_dir <- args[1]
out_dir   <- args[2]
tmp_base  <- file.path(out_dir, "tmp_scarHRD")

dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_base, showWarnings = FALSE, recursive = TRUE)

files  <- list.files(input_dir, pattern = "scarHRD_input.*\\.txt$", full.names = TRUE)
n      <- length(files)

cat("input_dir :", input_dir, "\n")
cat("out_dir   :", out_dir, "\n")
cat("n_files   :", n, "\n\n")

if (n == 0) {
  cat("ERROR: no scarHRD input files found\n")
  quit(status = 1)
}

original_wd    <- getwd()
all_results    <- list()
failed_samples <- character()

for (i in seq_along(files)) {
  f      <- files[i]
  sample <- sub("_scarHRD_input\\.txt$", "", basename(f))
  tmpdir <- file.path(tmp_base, sprintf("s%04d", i))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("[%d/%d] %s ... ", i, n, sample))

  tryCatch({
    setwd(tmpdir)
    f_tmp <- file.path(tmpdir, basename(f))
    file.copy(f, f_tmp, overwrite = TRUE)
    res <- scar_score(f_tmp, reference = "grch37", seqz = FALSE)
    df  <- as.data.frame(res)
    df$SampleID  <- sample
    df$InputFile <- basename(f)
    all_results[[sample]] <- df
    cat("OK\n")
  }, error = function(e) {
    cat("FAIL:", conditionMessage(e), "\n")
    failed_samples <<- c(failed_samples, sample)
  }, finally = {
    setwd(original_wd)
    unlink(tmpdir, recursive = TRUE, force = TRUE)
  })
}

unlink(tmp_base, recursive = TRUE, force = TRUE)

if (length(all_results) > 0) {
  final_df    <- do.call(rbind, all_results)
  output_file <- file.path(out_dir, "ALL_scarHRD_scores.csv")
  write.csv(final_df, file = output_file, row.names = FALSE)
  cat("\nProcessed:", length(all_results), "| Failed:", length(failed_samples), "\n")
  cat("Output:", output_file, "\n")
  print(final_df)
} else {
  cat("ERROR: no results generated\n")
  quit(status = 1)
}

if (length(failed_samples) > 0) {
  writeLines(failed_samples, file.path(out_dir, "failed_samples.txt"))
  cat("Failed samples written to failed_samples.txt\n")
}
