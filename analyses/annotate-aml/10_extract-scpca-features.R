#!/usr/bin/env Rscript
# Extract a compact per-cell table from the SCPCP000007 processed objects, for checking the
# BoneMarrowMap labels: ADT values for all 21 targets, log-normalized expression of a marker gene
# panel, QC columns, ScPCA's cell type annotations and ScPCA's per-cell inferCNV columns.
#
# See .workspace/annotate-aml/scpcp000007.md (step 3). The report
# (11_validate-scpca-labels.Rmd) reads this table, so it doesn't have to read 28 objects.
#
# Values are ScPCA's own normalization: `logcounts` for RNA and the ADT altExp's `logcounts`
# (median-based). ADT values are kept for every cell; `adt_scpca_filter` says whether ScPCA
# flagged the cell's ADT data.
#
# Input:
#   - data/scpca/SCPCP000007/single_cell_metadata.tsv
#   - data/scpca/SCPCP000007/<sample>/<library>_processed.rds
# Output (results/scpca/features/):
#   - libraries/<library>_features.tsv
#   - cell_features_SCPCP000007.tsv: all libraries combined
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 10_extract-scpca-features.R
# Libraries that already have output are skipped unless FORCE=1 is set.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

force <- Sys.getenv("FORCE", "0") == "1"

# Marker genes (RNA), by what they mark. Kept as columns `rna_<gene>`.
rna_markers <- c(
  # HSC/MPP
  "CD34", "CRHBP", "AVP", "HLF", "MECOM", "PROM1", "SPINK2",
  # lymphoid-primed progenitors
  "FLT3", "SATB1", "DNTT",
  # GMP / neutrophil-primed
  "MPO", "ELANE", "AZU1", "PRTN3", "CTSG", "CEBPA", "CEBPE",
  # monocyte
  "LYZ", "CD14", "FCGR3A", "S100A8", "S100A9", "CSF1R",
  # dendritic cells
  "CLEC10A", "FCER1A", "CD1C", "IRF8", "LILRA4", "IL3RA",
  # erythroid
  "GATA1", "KLF1", "TFRC", "CA1", "GYPA", "ALAS2", "HBB", "HBA1", "HBA2", "HBD",
  # megakaryocyte
  "ITGA2B", "PF4", "PPBP", "GP9",
  # eosinophil/basophil/mast
  "HDC", "CPA3", "MS4A2", "TPSAB1", "CLC", "GATA2",
  # B lineage
  "CD79A", "CD19", "MS4A1", "VPREB1", "IGLL1", "MZB1", "JCHAIN",
  # T/NK
  "CD3E", "CD4", "CD8A", "IL7R", "GNLY", "NKG7",
  # cycling, pan-leukocyte
  "MKI67", "TOP2A", "PTPRC"
)
stopifnot(!anyDuplicated(rna_markers))

# colData columns to keep, if present
coldata_keep <- c(
  "sum", "detected", "subsets_mito_percent", "prob_compromised", "scDblFinder_class",
  "scDblFinder_score", "cluster", "adt_scpca_filter", "altexps_adt_sum",
  "consensus_celltype_annotation", "singler_celltype_annotation",
  "cellassign_celltype_annotation", "scimilarity_celltype_annotation",
  "infercnv_total_cnv", "is_infercnv_reference"
)

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
project_dir <- fs::path(root_dir, "data", "scpca", "SCPCP000007")
metadata_file <- fs::path(project_dir, "single_cell_metadata.tsv")
output_dir <- fs::path(module_dir, "results", "scpca", "features")
library_dir <- fs::path(output_dir, "libraries")
fs::dir_create(library_dir)
combined_file <- fs::path(output_dir, "cell_features_SCPCP000007.tsv")

stopifnot("SCPCP000007 metadata not found" = fs::file_exists(metadata_file))

# Libraries in scope (same as 09) ---------------------------------------------------------------
libraries <- data.table::fread(metadata_file, data.table = FALSE) |>
  dplyr::filter(diagnosis %in% c("Acute myeloid leukemia", "Non-cancerous")) |>
  dplyr::transmute(
    sample_id = scpca_sample_id,
    library_id = scpca_library_id,
    file = fs::path(project_dir, scpca_sample_id, paste0(scpca_library_id, "_processed.rds")),
    output_file = fs::path(library_dir, paste0(scpca_library_id, "_features.tsv"))
  )
stopifnot(nrow(libraries) == 28)

to_run <- if (force) libraries else libraries[!fs::file_exists(libraries$output_file), ]
message(nrow(to_run), " libraries to extract")

# Extract --------------------------------------------------------------------------------------
for (i in seq_len(nrow(to_run))) {
  lib <- to_run[i, ]
  message("[", i, "/", nrow(to_run), "] ", lib$library_id)
  sce <- readRDS(lib$file)
  cell_ids <- paste(lib$library_id, colnames(sce), sep = "-")

  # RNA markers: first gene with each symbol
  symbols <- as.character(rowData(sce)$gene_symbol)
  idx <- match(rna_markers, symbols)
  if (anyNA(idx)) {
    warning(lib$library_id, ": marker genes not found: ",
            paste(rna_markers[is.na(idx)], collapse = ", "))
  }
  found <- !is.na(idx)
  rna <- as.matrix(logcounts(sce)[idx[found], , drop = FALSE])
  rna_df <- as.data.frame(t(rna))
  colnames(rna_df) <- paste0("rna_", rna_markers[found])

  # ADT
  stopifnot("adt" %in% altExpNames(sce))
  adt <- altExp(sce, "adt")
  adt_df <- as.data.frame(t(as.matrix(logcounts(adt))))
  colnames(adt_df) <- paste0("adt_", rownames(adt))

  cd <- as.data.frame(colData(sce))
  keep <- intersect(coldata_keep, colnames(cd))
  cd_df <- cd[, keep, drop = FALSE]
  for (col in setdiff(coldata_keep, keep)) cd_df[[col]] <- NA

  features <- data.frame(
    cell_id = cell_ids,
    library_id = lib$library_id,
    sample_id = lib$sample_id,
    barcode = colnames(sce),
    cd_df[, coldata_keep],
    rna_df,
    adt_df,
    check.names = FALSE,
    row.names = NULL
  )

  tmp_file <- paste0(lib$output_file, ".tmp")
  data.table::fwrite(features, tmp_file, sep = "\t", na = "NA")
  fs::file_move(tmp_file, lib$output_file)

  rm(sce, adt, rna, rna_df, adt_df, cd, cd_df, features)
  gc(verbose = FALSE)
}

# Combine --------------------------------------------------------------------------------------
combined <- purrr::map(libraries$output_file, \(f) data.table::fread(f, data.table = FALSE)) |>
  dplyr::bind_rows()
stopifnot(!anyDuplicated(combined$cell_id))
data.table::fwrite(combined, combined_file, sep = "\t", na = "NA")
message("Wrote ", combined_file, " (", nrow(combined), " cells, ", ncol(combined), " columns)")
