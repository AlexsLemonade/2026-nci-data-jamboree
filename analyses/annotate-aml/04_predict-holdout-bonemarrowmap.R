#!/usr/bin/env Rscript
# Annotate held-out atlas cells with BoneMarrowMap: Symphony mapping onto a normal bone marrow
# reference, mapping-error QC, then kNN cell type and pseudotime prediction.
#
# See .workspace/annotate-aml/bonemarrowmap.md (step 1).
#
# BoneMarrowMap's reference is built only from healthy bone marrow, so there's no held-out model
# to train; the held-out study is just the query.
#
# Input:
#   - data/AML-scAtlas/AML_scAtlas_le29.h5ad (01_subset-atlas.py)
#   - data/BoneMarrowMap/BoneMarrowMap_SymphonyReference.rds and
#     data/BoneMarrowMap/BoneMarrowMap_uwot_model.uwot (00_download-reference-data.sh)
# Output:
#   - results/holdout/bonemarrowmap_predictions_<study>.tsv: one row per cell, labels for every
#     cell (including mapping-QC failures; filter on bmm_mapping_qc downstream)
#   - results/holdout/bonemarrowmap_plots_<study>/: the package's per-donor projection plots
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 04_predict-holdout-bonemarrowmap.R [holdout_study]
# Default: holdout_study = mumme_2023

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BoneMarrowMap)
})

args <- commandArgs(trailingOnly = TRUE)
holdout_study <- if (length(args) >= 1) args[[1]] else "mumme_2023"

# Batch variable for Symphony's query correction (decision 2)
batch_var <- "Sample"
# Mapping-QC threshold, in MADs above the query median (decision 3: package default)
mad_threshold <- 2.5
# Stop if fewer than this fraction of the reference's variable genes are found in the query
min_gene_overlap <- 0.8

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
atlas_file <- fs::path(root_dir, "data", "AML-scAtlas", "AML_scAtlas_le29.h5ad")
bmm_dir <- fs::path(root_dir, "data", "BoneMarrowMap")
ref_file <- fs::path(bmm_dir, "BoneMarrowMap_SymphonyReference.rds")
uwot_file <- fs::path(bmm_dir, "BoneMarrowMap_uwot_model.uwot")
output_dir <- fs::path(module_dir, "results", "holdout")
fs::dir_create(output_dir)
output_file <- fs::path(output_dir, paste0("bonemarrowmap_predictions_", holdout_study, ".tsv"))
plot_dir <- fs::path(output_dir, paste0("bonemarrowmap_plots_", holdout_study))

stopifnot(
  "Atlas subset not found" = fs::file_exists(atlas_file),
  "BoneMarrowMap reference not found; run 00_download-reference-data.sh" =
    fs::file_exists(ref_file),
  "BoneMarrowMap uwot model not found; run 00_download-reference-data.sh" =
    fs::file_exists(uwot_file)
)

# Held-out cells -------------------------------------------------------------------------------
message("Reading ", holdout_study, " cells from ", atlas_file)
adata <- anndataR::read_h5ad(atlas_file, as = "HDF5AnnData")
is_holdout <- adata$obs$Study == holdout_study
stopifnot("Held-out study not found in the atlas" = any(is_holdout))
message(sum(is_holdout), " cells")

obs <- adata$obs[is_holdout, , drop = FALSE]
cell_ids <- adata$obs_names[is_holdout]
gene_symbols <- setNames(as.character(adata$var$gene_symbol), adata$var_names)

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
counts <- as(counts(test_sce), "CsparseMatrix")
rm(test_sce)
stopifnot(identical(colnames(counts), cell_ids))

# Atlas cells are already QC-filtered; log the ranges so that can be checked
n_count <- Matrix::colSums(counts)
n_feature <- Matrix::colSums(counts > 0)
message("nCount range: ", paste(range(n_count), collapse = "-"),
        "; nFeature range: ", paste(range(n_feature), collapse = "-"))

# Query metadata, indexed by cell ID (map_Query takes cell names from the row names)
metadata <- data.frame(
  Sample = as.character(obs$Sample),
  donor_id = as.character(obs$donor_id),
  assay = as.character(obs$assay),
  row.names = cell_ids
)

# Decision 2: the batch variable is Sample, unless a sample spans both assays
message("Cells by sample and assay:")
sample_assay <- table(metadata$Sample, metadata$assay)
print(sample_assay)
if (any(rowSums(sample_assay > 0) > 1)) {
  stop("Some samples have cells from more than one assay; decision 2 says to switch the ",
       "batch variable to Sample + assay. Update the plan and batch_var before rerunning.")
}

# Gene symbols ---------------------------------------------------------------------------------
# BoneMarrowMap matches genes by symbol. The atlas is indexed by Ensembl ID.
symbols <- unname(gene_symbols[rownames(counts)])
has_symbol <- !is.na(symbols) & symbols != ""
message("Dropping ", sum(!has_symbol), " of ", length(symbols), " genes with no symbol")
counts <- counts[has_symbol, ]
symbols <- symbols[has_symbol]

# Seurat replaces underscores in feature names with dashes, and map_Query uses the matrix's
# names in some places and the Seurat object's in others, so make them agree up front
has_underscore <- grepl("_", symbols, fixed = TRUE)
if (any(has_underscore)) {
  message("Replacing underscores with dashes in ", sum(has_underscore), " gene symbols")
  symbols <- gsub("_", "-", symbols, fixed = TRUE)
}

# For duplicated symbols, keep the gene with the most counts
gene_totals <- Matrix::rowSums(counts)
keep_order <- order(gene_totals, decreasing = TRUE)
keep <- keep_order[!duplicated(symbols[keep_order])]
keep <- sort(keep)
message("Dropping ", length(symbols) - length(keep), " genes with duplicated symbols")
counts <- counts[keep, ]
rownames(counts) <- symbols[keep]
message(nrow(counts), " genes in the query")

# Reference ------------------------------------------------------------------------------------
message("Reading ", ref_file)
ref <- readRDS(ref_file)

# The reference stores a relative uwot path with a different file name (decision 5)
ref$save_uwot_path <- as.character(uwot_file)

# map_Query reads ref$vargenes, which only resolves to ref$vargenes_means_sds through R's
# partial `$` matching on lists; check that it does
ref_genes <- ref$vargenes$symbol
stopifnot("Reference variable genes not found" = length(ref_genes) > 0)
if (any(grepl("_", ref_genes, fixed = TRUE))) {
  warning("Some reference variable genes contain underscores and won't match the query")
}
n_shared <- sum(ref_genes %in% rownames(counts))
message(n_shared, " of ", length(ref_genes), " reference variable genes found in the query (",
        round(100 * n_shared / length(ref_genes), 1), "%)")
if (n_shared / length(ref_genes) < min_gene_overlap) {
  stop("Gene overlap with the reference is below ", 100 * min_gene_overlap, "%")
}
message("Missing reference variable genes: ",
        paste(setdiff(ref_genes, rownames(counts)), collapse = ", "))

# Fine -> broad label lookup, used for every cell (the package only fills the broad label for
# cells that pass mapping QC)
label_lookup <- unique(ref$meta_data[c("CellType_Annotation", "CellType_Broad")])
stopifnot("A reference label has more than one broad label" =
            !anyDuplicated(label_lookup$CellType_Annotation))

# Map and predict ------------------------------------------------------------------------------
message("Mapping query cells (batch variable: ", batch_var, ")")
query <- map_Query(
  exp_query = counts,
  metadata_query = metadata,
  ref_obj = ref,
  vars = batch_var
)
rm(counts)

query <- calculate_MappingError(query, reference = ref, MAD_threshold = mad_threshold)
qc_counts <- table(query$mapping_error_QC, useNA = "ifany")
message("Mapping QC (", mad_threshold, " MADs, all cells): ",
        paste(names(qc_counts), qc_counts, sep = " = ", collapse = ", "))
message("Mapping QC failures by assay:")
print(table(query$assay, query$mapping_error_QC, useNA = "ifany"))

query <- predict_CellTypes(
  query_obj = query,
  ref_obj = ref,
  initial_label = "initial_CellType",
  final_label = "predicted_CellType"
)
query <- predict_Pseudotime(
  query_obj = query,
  ref_obj = ref,
  initial_label = "initial_Pseudotime",
  final_label = "predicted_Pseudotime"
)

# Write ----------------------------------------------------------------------------------------
meta <- query@meta.data
stopifnot(identical(rownames(meta), cell_ids))

if ("initial_CellType_prob" %in% colnames(meta)) {
  label_prob <- meta$initial_CellType_prob
} else {
  warning("initial_CellType_prob not found; bmm_label_prob will be NA")
  label_prob <- NA_real_
}

umap <- Seurat::Embeddings(query, reduction = "umap_projected")[cell_ids, , drop = FALSE]

# Labels and pseudotime are kept for every cell, including mapping-QC failures (decision 3).
# The package's QC-masked columns (predicted_*) are just these with failures set to NA.
pred_df <- data.frame(
  cell_id = cell_ids,
  donor_id = meta$donor_id,
  sample = meta$Sample,
  assay = meta$assay,
  bmm_label = meta$initial_CellType,
  bmm_label_broad = label_lookup$CellType_Broad[
    match(meta$initial_CellType, label_lookup$CellType_Annotation)
  ],
  bmm_label_prob = label_prob,
  bmm_pseudotime = meta$initial_Pseudotime,
  bmm_mapping_error = meta$mapping_error_score,
  bmm_mapping_qc = meta$mapping_error_QC,
  bmm_umap_1 = umap[, 1],
  bmm_umap_2 = umap[, 2]
)
stopifnot(!anyNA(pred_df$bmm_label), !anyNA(pred_df$bmm_label_broad))
data.table::fwrite(pred_df, output_file, sep = "\t", na = "NA")
message("Wrote ", output_file)

message("Cells by BoneMarrowMap broad label:")
print(sort(table(pred_df$bmm_label_broad), decreasing = TRUE))

# Per-donor projection plots -------------------------------------------------------------------
# Not needed for the TSV, so a failure here (e.g. from Seurat 6) skips the plots with a message
fs::dir_create(plot_dir)
tryCatch(
  {
    plot_Projection_byDonor(
      query_obj = query,
      batch_key = "donor_id",
      ref_obj = ref,
      save_folder = paste0(plot_dir, "/")
    )
    message("Wrote per-donor projection plots to ", plot_dir)
  },
  error = function(e) {
    message("Skipping per-donor projection plots: ", conditionMessage(e))
  }
)
