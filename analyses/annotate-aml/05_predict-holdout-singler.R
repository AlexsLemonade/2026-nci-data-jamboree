#!/usr/bin/env Rscript
# Classify held-out atlas cells with the SingleR reference trained without them.
#
# See .workspace/annotate-aml/evaluation.md (design step 1).
#
# Input:
#   - data/AML-scAtlas/AML_scAtlas_le29.h5ad (02_subset-atlas.py)
#   - results/AML_scAtlas_le29_singler_trained_ref_no-<study>.rds (03_build-singler-ref.Rmd
#     rendered with holdout_study = <study>)
# Output:
#   - results/holdout/singler_predictions_<study>.tsv
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 05_predict-holdout-singler.R [holdout_study] [n_workers]
# Defaults: holdout_study = mumme_2023, n_workers = 4

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SingleR)
})

args <- commandArgs(trailingOnly = TRUE)
holdout_study <- if (length(args) >= 1) args[[1]] else "mumme_2023"
n_workers <- if (length(args) >= 2) as.integer(args[[2]]) else 4L

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
atlas_file <- fs::path(root_dir, "data", "AML-scAtlas", "AML_scAtlas_le29.h5ad")
trained_ref_file <- fs::path(
  module_dir, "results",
  paste0("AML_scAtlas_le29_singler_trained_ref_no-", holdout_study, ".rds")
)
output_dir <- fs::path(module_dir, "results", "holdout")
fs::dir_create(output_dir)
output_file <- fs::path(output_dir, paste0("singler_predictions_", holdout_study, ".tsv"))

stopifnot(
  "Atlas subset not found" = fs::file_exists(atlas_file),
  "Held-out trained reference not found; render 03_build-singler-ref.Rmd with holdout_study" =
    fs::file_exists(trained_ref_file)
)

# Held-out cells -------------------------------------------------------------------------------
message("Reading ", holdout_study, " cells from ", atlas_file)
adata <- anndataR::read_h5ad(atlas_file, as = "HDF5AnnData")
is_holdout <- adata$obs$Study == holdout_study
stopifnot("Held-out study not found in the atlas" = any(is_holdout))
message(sum(is_holdout), " cells")

test_sce <- adata[is_holdout, ]$as_SingleCellExperiment(
  x_mapping = "counts",
  assays_mapping = FALSE,
  colData_mapping = FALSE,
  rowData_mapping = FALSE,
  reducedDims_mapping = FALSE,
  colPairs_mapping = FALSE,
  rowPairs_mapping = FALSE,
  metadata_mapping = FALSE
)
rm(adata)
counts(test_sce) <- as(counts(test_sce), "CsparseMatrix")

# Same normalization as the reference
test_sce <- scuttle::logNormCounts(test_sce)

# Classify -------------------------------------------------------------------------------------
trained_ref <- readRDS(trained_ref_file)

message("Classifying with ", trained_ref_file)
pred <- classifySingleR(
  test = test_sce,
  trained = trained_ref,
  BPPARAM = BiocParallel::MulticoreParam(workers = n_workers)
)

message(
  sum(is.na(pred$pruned.labels)), " of ", nrow(pred), " cells have pruned (NA) labels"
)

# Write ----------------------------------------------------------------------------------------
pred_df <- data.frame(
  cell_id = colnames(test_sce),
  singler_label = pred$labels,
  singler_pruned_label = pred$pruned.labels,
  singler_delta_next = pred$delta.next
)
data.table::fwrite(pred_df, output_file, sep = "\t", na = "NA")
message("Wrote ", output_file)
