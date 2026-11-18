#!/usr/bin/env Rscript
# Annotate SCPCP000007 cells with BoneMarrowMap, one library at a time: Symphony mapping onto
# the normal bone marrow reference, mapping-error scores, kNN cell type and pseudotime.
#
# See .workspace/annotate-aml/scpcp000007.md (step 2). Follows 04_predict-holdout-bonemarrowmap.R.
#
# Libraries: the 26 AML libraries and the 2 healthy controls (the healthy controls are the
# normal reference for later steps). Each library is one sample, so no batch variable is used.
#
# Mapping QC: calculate_MappingError() flags cells relative to the query it's given, which here
# is one library. That per-library flag is saved as `bmm_mapping_qc_library`, but the plan uses a
# threshold across the whole project (decision 3 in bonemarrowmap.md), computed later from
# `bmm_mapping_error`.
#
# Input:
#   - data/scpca/SCPCP000007/single_cell_metadata.tsv
#   - data/scpca/SCPCP000007/<sample>/<library>_processed.rds
#   - data/BoneMarrowMap/BoneMarrowMap_SymphonyReference.rds and
#     data/BoneMarrowMap/BoneMarrowMap_uwot_model.uwot (00_download-reference-data.sh)
# Output (results/scpca/bonemarrowmap/):
#   - libraries/<library>_bonemarrowmap.tsv: one row per cell, written as each library finishes
#   - library_summary.tsv: cells, gene overlap and per-library QC failures for each library
#   - bonemarrowmap_predictions_SCPCP000007.tsv: all libraries combined
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 09_predict-scpca-bonemarrowmap.R [library_id ...]
# With no arguments, all libraries in scope are run. Libraries that already have output are
# skipped unless FORCE=1 is set.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BoneMarrowMap)
})

args <- commandArgs(trailingOnly = TRUE)
force <- Sys.getenv("FORCE", "0") == "1"

# Same settings as 04
mad_threshold <- 2.5
min_gene_overlap <- 0.8

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
project_dir <- fs::path(root_dir, "data", "scpca", "SCPCP000007")
metadata_file <- fs::path(project_dir, "single_cell_metadata.tsv")
bmm_dir <- fs::path(root_dir, "data", "BoneMarrowMap")
ref_file <- fs::path(bmm_dir, "BoneMarrowMap_SymphonyReference.rds")
uwot_file <- fs::path(bmm_dir, "BoneMarrowMap_uwot_model.uwot")
output_dir <- fs::path(module_dir, "results", "scpca", "bonemarrowmap")
library_dir <- fs::path(output_dir, "libraries")
fs::dir_create(library_dir)
summary_file <- fs::path(output_dir, "library_summary.tsv")
combined_file <- fs::path(output_dir, "bonemarrowmap_predictions_SCPCP000007.tsv")

stopifnot(
  "SCPCP000007 metadata not found" = fs::file_exists(metadata_file),
  "BoneMarrowMap reference not found; run 00_download-reference-data.sh" =
    fs::file_exists(ref_file),
  "BoneMarrowMap uwot model not found; run 00_download-reference-data.sh" =
    fs::file_exists(uwot_file)
)

# Libraries in scope ---------------------------------------------------------------------------
libraries <- data.table::fread(metadata_file, data.table = FALSE) |>
  dplyr::filter(diagnosis %in% c("Acute myeloid leukemia", "Non-cancerous")) |>
  dplyr::transmute(
    sample_id = scpca_sample_id,
    library_id = scpca_library_id,
    role = dplyr::if_else(diagnosis == "Non-cancerous", "healthy control", "AML"),
    file = fs::path(project_dir, scpca_sample_id, paste0(scpca_library_id, "_processed.rds")),
    output_file = fs::path(library_dir, paste0(scpca_library_id, "_bonemarrowmap.tsv"))
  )
stopifnot(sum(libraries$role == "AML") == 26, sum(libraries$role == "healthy control") == 2)

if (length(args) > 0) {
  unknown <- setdiff(args, libraries$library_id)
  if (length(unknown) > 0) stop("Libraries not in scope: ", paste(unknown, collapse = ", "))
  to_run <- libraries[libraries$library_id %in% args, ]
} else {
  to_run <- libraries
}
if (!force) {
  done <- fs::file_exists(to_run$output_file)
  if (any(done)) {
    message("Skipping ", sum(done), " libraries with existing output (set FORCE=1 to rerun)")
  }
  to_run <- to_run[!done, ]
}
message(nrow(to_run), " libraries to map")

# Reference ------------------------------------------------------------------------------------
message("Reading ", ref_file)
ref <- readRDS(ref_file)
# The reference stores a relative uwot path with a different file name (decision 5 in
# bonemarrowmap.md)
ref$save_uwot_path <- as.character(uwot_file)
# map_Query reads ref$vargenes, which resolves to ref$vargenes_means_sds by partial matching
ref_genes <- ref$vargenes$symbol
stopifnot("Reference variable genes not found" = length(ref_genes) > 0)

label_lookup <- unique(ref$meta_data[c("CellType_Annotation", "CellType_Broad")])
stopifnot("A reference label has more than one broad label" =
            !anyDuplicated(label_lookup$CellType_Annotation))

# Counts with gene symbols as row names, prepared as in 04: drop genes with no symbol, replace
# underscores with dashes, and for duplicated symbols keep the gene with the most counts
prepare_counts <- function(sce) {
  counts <- as(counts(sce), "CsparseMatrix")
  symbols <- as.character(rowData(sce)$gene_symbol)
  has_symbol <- !is.na(symbols) & symbols != ""
  counts <- counts[has_symbol, ]
  symbols <- gsub("_", "-", symbols[has_symbol], fixed = TRUE)
  gene_totals <- Matrix::rowSums(counts)
  keep_order <- order(gene_totals, decreasing = TRUE)
  keep <- sort(keep_order[!duplicated(symbols[keep_order])])
  counts <- counts[keep, ]
  rownames(counts) <- symbols[keep]
  list(counts = counts, n_no_symbol = sum(!has_symbol), n_duplicated = length(symbols) - length(keep))
}

# Map each library -----------------------------------------------------------------------------
for (i in seq_len(nrow(to_run))) {
  lib <- to_run[i, ]
  message("\n[", i, "/", nrow(to_run), "] ", lib$library_id, " (", lib$sample_id, ", ",
          lib$role, ")")
  sce <- readRDS(lib$file)
  barcodes <- colnames(sce)
  stopifnot("Barcodes not unique" = !anyDuplicated(barcodes))
  cell_ids <- paste(lib$library_id, barcodes, sep = "-")

  prepared <- prepare_counts(sce)
  rm(sce)
  counts <- prepared$counts
  colnames(counts) <- cell_ids

  n_shared <- sum(ref_genes %in% rownames(counts))
  message(ncol(counts), " cells; ", n_shared, " of ", length(ref_genes),
          " reference variable genes found (", round(100 * n_shared / length(ref_genes), 1), "%)")
  if (n_shared / length(ref_genes) < min_gene_overlap) {
    stop("Gene overlap with the reference is below ", 100 * min_gene_overlap, "%")
  }

  metadata <- data.frame(library_id = rep(lib$library_id, length(cell_ids)), row.names = cell_ids)

  query <- map_Query(exp_query = counts, metadata_query = metadata, ref_obj = ref, vars = NULL,
                     verbose = FALSE)
  rm(counts)
  query <- calculate_MappingError(query, reference = ref, MAD_threshold = mad_threshold)
  query <- predict_CellTypes(query_obj = query, ref_obj = ref, initial_label = "initial_CellType",
                             final_label = "predicted_CellType")
  query <- predict_Pseudotime(query_obj = query, ref_obj = ref,
                              initial_label = "initial_Pseudotime",
                              final_label = "predicted_Pseudotime")

  meta <- query@meta.data
  stopifnot(identical(rownames(meta), cell_ids))
  umap <- Seurat::Embeddings(query, reduction = "umap_projected")[cell_ids, , drop = FALSE]
  rm(query)

  label_prob <- if ("initial_CellType_prob" %in% colnames(meta)) {
    meta$initial_CellType_prob
  } else {
    warning("initial_CellType_prob not found; bmm_label_prob will be NA")
    NA_real_
  }

  # As in 04: labels and pseudotime for every cell (the package's initial_* columns); the raw
  # pseudotime is kept, and off-hierarchy labels are handled downstream
  pred_df <- data.frame(
    cell_id = cell_ids,
    library_id = lib$library_id,
    sample_id = lib$sample_id,
    barcode = barcodes,
    role = lib$role,
    bmm_label = meta$initial_CellType,
    bmm_label_broad = label_lookup$CellType_Broad[
      match(meta$initial_CellType, label_lookup$CellType_Annotation)
    ],
    bmm_label_prob = label_prob,
    bmm_pseudotime = meta$initial_Pseudotime,
    bmm_mapping_error = meta$mapping_error_score,
    bmm_mapping_qc_library = meta$mapping_error_QC,
    bmm_umap_1 = umap[, 1],
    bmm_umap_2 = umap[, 2]
  )
  stopifnot(!anyNA(pred_df$bmm_label), !anyNA(pred_df$bmm_label_broad))

  # The summary is written first and the predictions last (through a temporary file), so a
  # library only counts as done once both exist
  data.table::fwrite(
    data.frame(
      library_id = lib$library_id, sample_id = lib$sample_id, role = lib$role,
      n_cells = nrow(pred_df), n_genes_no_symbol = prepared$n_no_symbol,
      n_genes_duplicated_symbol = prepared$n_duplicated, n_ref_genes_found = n_shared,
      n_ref_genes = length(ref_genes),
      n_fail_library_qc = sum(pred_df$bmm_mapping_qc_library == "Fail")
    ),
    fs::path(library_dir, paste0(lib$library_id, "_summary.tsv")),
    sep = "\t"
  )
  tmp_file <- paste0(lib$output_file, ".tmp")
  data.table::fwrite(pred_df, tmp_file, sep = "\t", na = "NA")
  fs::file_move(tmp_file, lib$output_file)
  message("Wrote ", lib$output_file, " (", sum(pred_df$bmm_mapping_qc_library == "Fail"),
          " cells fail per-library QC)")

  rm(pred_df, meta, umap, prepared)
  gc(verbose = FALSE)
}

# Combine --------------------------------------------------------------------------------------
# Only once every library in scope has output
all_done <- all(fs::file_exists(libraries$output_file))
if (all_done) {
  combined <- purrr::map(libraries$output_file, \(f) data.table::fread(f, data.table = FALSE)) |>
    dplyr::bind_rows()
  stopifnot(!anyDuplicated(combined$cell_id))
  data.table::fwrite(combined, combined_file, sep = "\t", na = "NA")
  message("\nWrote ", combined_file, " (", nrow(combined), " cells, ",
          dplyr::n_distinct(combined$library_id), " libraries)")

  summary_files <- fs::path(library_dir, paste0(libraries$library_id, "_summary.tsv"))
  library_summary <- purrr::map(summary_files, \(f) data.table::fread(f, data.table = FALSE)) |>
    dplyr::bind_rows()
  data.table::fwrite(library_summary, summary_file, sep = "\t")
  message("Wrote ", summary_file)
  print(library_summary, row.names = FALSE)

  message("\nCells by BoneMarrowMap broad label and role:")
  print(table(combined$bmm_label_broad, combined$role))
} else {
  message("\nNot all libraries have output yet; skipping the combined files")
}
