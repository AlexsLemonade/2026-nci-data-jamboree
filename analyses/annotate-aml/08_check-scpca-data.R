#!/usr/bin/env Rscript
# Read-only check of the SCPCP000007 processed objects used for annotation: the AML libraries
# and the healthy controls (the T-myeloid MPAL samples are out of scope).
#
# See .workspace/annotate-aml/scpcp000007.md (step 1).
#
# For each library's processed SCE (.rds), records:
#   - cell and gene counts, assay names, whether `counts` holds raw (integer) counts
#   - whether rowData has gene symbols
#   - the ADT altExp: name, targets, assays
#   - ScPCA's inferCNV status and reference (from metadata) and the per-cell inferCNV columns
#   - the consensus cell type labels
#
# Input:
#   - data/scpca/SCPCP000007/single_cell_metadata.tsv
#   - data/scpca/SCPCP000007/<sample>/<library>_processed.rds
# Output (results/scpca/data-check/):
#   - library_summary.tsv: one row per library
#   - coldata_columns.tsv: which colData columns each library has
#   - adt_targets.tsv: ADT targets per library
#   - consensus_labels.tsv: consensus cell type counts per library
#   - infercnv_reference.tsv: ScPCA's inferCNV reference cell types per library
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 08_check-scpca-data.R

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
project_dir <- fs::path(root_dir, "data", "scpca", "SCPCP000007")
metadata_file <- fs::path(project_dir, "single_cell_metadata.tsv")
output_dir <- fs::path(module_dir, "results", "scpca", "data-check")
fs::dir_create(output_dir)

stopifnot("SCPCP000007 metadata not found" = fs::file_exists(metadata_file))

# Libraries in scope ---------------------------------------------------------------------------
# AML libraries are labeled; healthy controls are the normal reference (scpcp000007.md, Scope)
sample_metadata <- data.table::fread(metadata_file, data.table = FALSE)
libraries <- sample_metadata |>
  dplyr::filter(diagnosis %in% c("Acute myeloid leukemia", "Non-cancerous")) |>
  dplyr::transmute(
    sample_id = scpca_sample_id,
    library_id = scpca_library_id,
    role = dplyr::if_else(diagnosis == "Non-cancerous", "healthy control", "AML"),
    diagnosis,
    subdiagnosis,
    disease_timing,
    processed_cells_metadata = processed_cells,
    file = fs::path(project_dir, scpca_sample_id, paste0(scpca_library_id, "_processed.rds"))
  )
message(nrow(libraries), " libraries: ", paste(names(table(libraries$role)), table(libraries$role),
                                              sep = " = ", collapse = ", "))
stopifnot("Expected 26 AML libraries and 2 healthy controls" =
            sum(libraries$role == "AML") == 26 && sum(libraries$role == "healthy control") == 2)

missing <- libraries$file[!fs::file_exists(libraries$file)]
if (length(missing) > 0) stop("Missing processed files:\n", paste(missing, collapse = "\n"))

# Helpers --------------------------------------------------------------------------------------
# Collapse a metadata value (possibly a vector or NULL) to one string
as_string <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  paste(as.character(unlist(x)), collapse = "; ")
}

# Whether a sparse or dense matrix holds whole numbers (checks up to `n` non-zero values)
is_integer_valued <- function(mat, n = 1e6) {
  values <- if (methods::is(mat, "dgCMatrix")) mat@x else as.vector(as.matrix(mat))
  values <- values[values != 0]
  values <- head(values, n)
  length(values) > 0 && all(values == round(values))
}

infercnv_cell_columns <- c("infercnv_total_cnv", "is_infercnv_reference")

# Check each library ---------------------------------------------------------------------------
summaries <- list()
coldata_columns <- list()
adt_targets <- list()
consensus_labels <- list()
infercnv_reference <- list()

for (i in seq_len(nrow(libraries))) {
  lib <- libraries[i, ]
  message("\n[", i, "/", nrow(libraries), "] ", lib$library_id, " (", lib$sample_id, ", ",
          lib$role, ")")
  sce <- readRDS(lib$file)

  # Print the full metadata names once, to see what else ScPCA stores
  if (i == 1) {
    message("metadata names: ", paste(names(metadata(sce)), collapse = ", "))
    message("rowData columns: ", paste(colnames(rowData(sce)), collapse = ", "))
    message("assays: ", paste(assayNames(sce), collapse = ", "),
            "; altExps: ", paste(altExpNames(sce), collapse = ", "))
  }

  cd <- colData(sce)
  meta <- metadata(sce)

  # ADT
  adt_name <- intersect(c("adt", "ADT"), altExpNames(sce))
  adt_name <- if (length(adt_name) > 0) adt_name[[1]] else NA_character_
  adt <- if (!is.na(adt_name)) altExp(sce, adt_name) else NULL
  if (!is.null(adt)) {
    adt_targets[[lib$library_id]] <- data.frame(library_id = lib$library_id, target = rownames(adt))
  }

  # inferCNV metadata: any entry whose name mentions infercnv
  infercnv_meta <- meta[grepl("infercnv", names(meta), ignore.case = TRUE)]
  if (i == 1 || !"infercnv_status" %in% names(infercnv_meta)) {
    message("inferCNV metadata entries: ", paste(names(infercnv_meta), collapse = ", "))
  }
  ref_celltypes <- meta[["infercnv_reference_celltypes"]]
  if (!is.null(ref_celltypes) && length(ref_celltypes) > 0) {
    infercnv_reference[[lib$library_id]] <- data.frame(
      library_id = lib$library_id, reference_celltype = as.character(unlist(ref_celltypes))
    )
  }

  # Consensus labels
  if ("consensus_celltype_annotation" %in% colnames(cd)) {
    tab <- table(cd$consensus_celltype_annotation, useNA = "ifany")
    consensus_labels[[lib$library_id]] <- data.frame(
      library_id = lib$library_id, consensus_celltype_annotation = names(tab),
      n_cells = as.integer(tab)
    )
  }

  coldata_columns[[lib$library_id]] <- data.frame(library_id = lib$library_id, column = colnames(cd))

  has_total_cnv <- "infercnv_total_cnv" %in% colnames(cd)
  has_ref_flag <- "is_infercnv_reference" %in% colnames(cd)

  summaries[[lib$library_id]] <- data.frame(
    sample_id = lib$sample_id,
    library_id = lib$library_id,
    role = lib$role,
    diagnosis = lib$diagnosis,
    subdiagnosis = lib$subdiagnosis,
    disease_timing = lib$disease_timing,
    n_cells = ncol(sce),
    n_cells_matches_metadata = ncol(sce) == lib$processed_cells_metadata,
    n_genes = nrow(sce),
    assays = paste(assayNames(sce), collapse = ","),
    counts_integer = "counts" %in% assayNames(sce) && is_integer_valued(counts(sce)),
    has_gene_symbol = "gene_symbol" %in% colnames(rowData(sce)),
    n_missing_gene_symbol = if ("gene_symbol" %in% colnames(rowData(sce)))
      sum(is.na(rowData(sce)$gene_symbol) | rowData(sce)$gene_symbol == "") else NA_integer_,
    adt_altexp = adt_name,
    n_adt_targets = if (!is.null(adt)) nrow(adt) else NA_integer_,
    adt_assays = if (!is.null(adt)) paste(assayNames(adt), collapse = ",") else NA_character_,
    n_consensus_labels = if ("consensus_celltype_annotation" %in% colnames(cd))
      length(unique(cd$consensus_celltype_annotation)) else NA_integer_,
    infercnv_status = as_string(meta[["infercnv_status"]]),
    infercnv_num_reference_cells = as_string(meta[["infercnv_num_reference_cells"]]),
    infercnv_diagnosis_groups = as_string(meta[["infercnv_diagnosis_groups"]]),
    has_infercnv_total_cnv = has_total_cnv,
    n_infercnv_reference_cells = if (has_ref_flag) sum(cd$is_infercnv_reference %in% TRUE) else
      NA_integer_,
    infercnv_total_cnv_median = if (has_total_cnv) median(cd$infercnv_total_cnv, na.rm = TRUE) else
      NA_real_,
    infercnv_total_cnv_frac_gt0 = if (has_total_cnv) mean(cd$infercnv_total_cnv > 0, na.rm = TRUE) else
      NA_real_,
    n_has_cnv_arm_columns = sum(grepl("^has_cnv_", colnames(cd)))
  )

  rm(sce, adt, cd, meta)
  gc(verbose = FALSE)
}

# Write ----------------------------------------------------------------------------------------
library_summary <- dplyr::bind_rows(summaries)
coldata_columns <- dplyr::bind_rows(coldata_columns)
adt_targets <- dplyr::bind_rows(adt_targets)
consensus_labels <- dplyr::bind_rows(consensus_labels)
infercnv_reference <- dplyr::bind_rows(infercnv_reference)

write_tsv <- function(df, name) {
  file <- fs::path(output_dir, name)
  data.table::fwrite(df, file, sep = "\t", na = "NA")
  message("Wrote ", file)
}
write_tsv(library_summary, "library_summary.tsv")
write_tsv(coldata_columns, "coldata_columns.tsv")
write_tsv(adt_targets, "adt_targets.tsv")
write_tsv(consensus_labels, "consensus_labels.tsv")
write_tsv(infercnv_reference, "infercnv_reference.tsv")

# Summary to the log ---------------------------------------------------------------------------
message("\nTotal cells: ", sum(library_summary$n_cells), " (AML: ",
        sum(library_summary$n_cells[library_summary$role == "AML"]), ")")
message("All counts integer-valued: ", all(library_summary$counts_integer))
message("All libraries have gene symbols: ", all(library_summary$has_gene_symbol))

# Is the ADT panel the same in every library?
adt_sets <- split(adt_targets$target, adt_targets$library_id)
same_adt <- length(unique(lapply(adt_sets, sort))) == 1
message("ADT panel identical across libraries: ", same_adt, " (",
        paste(sort(unique(adt_targets$target)), collapse = ", "), ")")

# colData columns not shared by every library
column_counts <- table(coldata_columns$column)
partial <- names(column_counts)[column_counts < nrow(library_summary)]
message("colData columns missing from some libraries: ",
        if (length(partial) > 0) paste(partial, collapse = ", ") else "none")

message("\nPer library:")
print(library_summary[, c("library_id", "role", "n_cells", "counts_integer", "n_adt_targets",
                          "infercnv_status", "infercnv_num_reference_cells",
                          "infercnv_total_cnv_frac_gt0")], row.names = FALSE)

message("\nConsensus labels, all libraries (top 30):")
consensus_labels |>
  dplyr::group_by(consensus_celltype_annotation) |>
  dplyr::summarize(n_cells = sum(n_cells), n_libraries = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_cells)) |>
  head(30) |>
  print()

message("\ninferCNV reference cell types:")
if (nrow(infercnv_reference) > 0) print(table(infercnv_reference$reference_celltype)) else
  message("none recorded")
