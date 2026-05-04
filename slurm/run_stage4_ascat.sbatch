#!/bin/bash
#SBATCH --job-name=scarHRD_ascat
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jmroldan@cnio.es

# ============================================================
# run_stage4_ascat.sbatch
#
# scarHRD scoring for the 20260406 ASCAT-downsampled samples.
# Reads Stage 3 outputs from scarHRD_inputs_thresholds/ and
# writes ALL_scarHRD_scores.csv to outputs_scarHRD_thresholds/.
#
# Submit after run_stage3_ascat.sbatch completes:
#   sbatch --dependency=afterok:<stage3_jobid> run_stage4_ascat.sbatch
# ============================================================

OUT_DIR="/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/outputs_scarHRD_thresholds"
mkdir -p "$OUT_DIR"

echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Output directory: $OUT_DIR"
echo "========================================"

conda run -n r43 R --vanilla --quiet <<'EOF'

library(scarHRD)

input_scarFiles <- "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/scarHRD_inputs_thresholds"
out_dir         <- "/storage/scratch01/groups/co/cn_extra/alleleSpecific/tso_samples/20260406/outputs_scarHRD_thresholds"
tmp_base_dir    <- file.path(out_dir, "tmp_processing")

dir.create(out_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_base_dir, showWarnings = FALSE, recursive = TRUE)

original_wd <- getwd()

# Each input file now contains all DP filters × both merge states for one
# sample×backbone, with columns: SampleID Chromosome Start_position End_position
# total_cn A_cn B_cn ploidy dp_filter merge_state.
# We iterate over files, then over (dp_filter, merge_state) combinations within
# each file, writing a slim temp file for scar_score().

files_scar <- list.files(
  input_scarFiles,
  pattern = "_scarHRD_input\\.txt$",
  full.names = TRUE,
  recursive = TRUE
)
n_files <- length(files_scar)

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║        scarHRD — ASCAT-downsampled 20260406                ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("Found", n_files, "input files\n\n")

if (n_files == 0) {
  cat("ERROR: No input files found!\n")
  cat("   Check path:", input_scarFiles, "\n")
  quit(status = 1)
}

scar_cols <- c("SampleID", "Chromosome", "Start_position", "End_position",
               "total_cn", "A_cn", "B_cn", "ploidy")

all_results    <- list()
failed_samples <- character()
start_time     <- Sys.time()
combo_idx      <- 0L

for (i in seq_along(files_scar)) {
  f     <- files_scar[i]
  fname <- basename(f)
  base_sample <- sub("_scarHRD_input\\.txt$", "", fname)

  full_df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) { cat(sprintf("  ✗ Cannot read %s: %s\n", fname, conditionMessage(e))); NULL }
  )
  if (is.null(full_df)) next

  dp_values    <- sort(unique(full_df$dp_filter))
  merge_states <- unique(full_df$merge_state)

  for (dp in dp_values) {
    for (ms in merge_states) {
      combo_idx <- combo_idx + 1L
      sample_id <- sprintf("%s_dp%s_%s", base_sample, dp, ms)

      sub_df <- full_df[full_df$dp_filter == dp & full_df$merge_state == ms, scar_cols, drop = FALSE]
      if (nrow(sub_df) == 0) next

      sample_tmp_dir <- file.path(tmp_base_dir, sprintf("combo_%08d", combo_idx))
      dir.create(sample_tmp_dir, showWarnings = FALSE, recursive = TRUE)

      cat(sprintf("[file %d/%d | dp=%s | %s] %s", i, n_files, dp, ms, base_sample))

      tryCatch({
        setwd(sample_tmp_dir)
        f_tmp <- file.path(sample_tmp_dir, "scarHRD_input.txt")
        write.table(sub_df, f_tmp, sep = "\t", row.names = FALSE, quote = FALSE)
        res <- scar_score(f_tmp, reference = "grch37", seqz = FALSE)
        df_res <- as.data.frame(res)
        df_res$SampleID    <- sample_id
        df_res$InputFile   <- f
        df_res$dp_filter   <- dp
        df_res$merge_state <- ms
        all_results[[sample_id]] <- df_res
        cat(" ✓")
      }, error = function(e) {
        cat(" ✗ Error:", conditionMessage(e))
        failed_samples <<- c(failed_samples, sample_id)
      }, finally = {
        setwd(original_wd)
      })

      gc()
      unlink(sample_tmp_dir, recursive = TRUE, force = TRUE)
      cat("\n")
    }
  }

  if (i %% 50 == 0) {
    elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate      <- combo_idx / elapsed
    n_combos  <- n_files * length(dp_values) * length(merge_states)
    remaining <- (n_combos - combo_idx) / rate
    stray_files <- list.files(out_dir, pattern = "_HRD.*\\.txt$|_chromosome.*\\.jpg$", full.names = FALSE)
    stray_files <- stray_files[!grepl("^scarHRD_[0-9]+\\.(log|err)$|ALL_scarHRD_scores\\.csv|failed_samples\\.txt", stray_files)]
    cat(sprintf("Progress: %d/%d files | ETA: %.1f min | Stray: %d\n",
                i, n_files, remaining, length(stray_files)))
    if (length(stray_files) > 0) file.remove(file.path(out_dir, stray_files))
  }

  if (i %% 200 == 0) {
    gc()
    unlink(tmp_base_dir, recursive = TRUE, force = TRUE)
    dir.create(tmp_base_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

setwd(original_wd)
end_time   <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

if (length(all_results) > 0) {
  final_df    <- do.call(rbind, all_results)
  output_file <- file.path(out_dir, "ALL_scarHRD_scores.csv")
  write.csv(final_df, file = output_file, row.names = FALSE)
  cat("\n✓ Processed:", length(all_results), "samples\n")
  cat("✓ Failed:",    length(failed_samples), "samples\n")
  cat("✓ Time:",      sprintf("%.1f", total_time), "minutes\n")
  cat("✓ Output:",    output_file, "\n")
} else {
  cat("\nERROR: No results generated\n")
}

if (length(failed_samples) > 0) {
  writeLines(failed_samples, file.path(out_dir, "failed_samples.txt"))
}

gc()
unlink(tmp_base_dir, recursive = TRUE, force = TRUE)
stray_files <- list.files(out_dir, pattern = "_HRD|_chromosome.*\\.jpg$|\\.pdf$", full.names = TRUE)
stray_files <- stray_files[!grepl("ALL_scarHRD_scores\\.csv|failed_samples\\.txt|scarHRD_[0-9]+\\.(log|err)$", basename(stray_files))]
if (length(stray_files) > 0) file.remove(stray_files)

cat("Done!\n")
EOF

echo ""
echo "========================================"
echo "End time: $(date)"
echo "========================================"
