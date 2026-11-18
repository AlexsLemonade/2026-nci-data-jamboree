#!/usr/bin/env Rscript
# Extract ScPCA's per-cell inferCNV results, including the arm-level calls, for the SCPCP000007
# AML libraries where ScPCA's inferCNV succeeded.
#
# See .workspace/annotate-aml/scpcp000007.md (step 4; decision 5). In the processed AnnData:
#   - `obs` has only `infercnv_total_cnv` and `is_infercnv_reference`;
#   - the arm-level HMM results are a separate data frame, `uns/infercnv_table`, keyed by
#     `barcodes`: per arm `has_cnv_*`, `has_loss_*`, `has_dupli_*`, `proportion_{cnv,loss,dupli}_*`,
#     `proportion_scaled_*`, plus `subcluster`.
# Neither is in the .rds. `uns/infercnv_table` is read with rhdf5 (anndataR's support for data
# frames in `uns` isn't relied on). Kept: `subcluster`, the `has_*` flags and the
# `proportion_loss_*` / `proportion_dupli_*` columns. Cells are matched to the other tables by
# `<library_id>-<barcode>`.
#
# Input:
#   - results/scpca/data-check/library_summary.tsv (08_check-scpca-data.R): inferCNV status
#   - data/scpca/SCPCP000007/<sample>/<library>_processed_rna.h5ad
#   - results/scpca/features/cell_features_SCPCP000007.tsv (10_extract-scpca-features.R): to check
#     that the barcodes match
# Output (results/scpca/infercnv/):
#   - scpca_infercnv_cells_SCPCP000007.tsv: one row per cell
#   - column_summary.tsv: type and value summary of each column, per library
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 12_extract-scpca-infercnv.R

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
project_dir <- fs::path(root_dir, "data", "scpca", "SCPCP000007")
scpca_results <- fs::path(module_dir, "results", "scpca")
library_summary_file <- fs::path(scpca_results, "data-check", "library_summary.tsv")
features_file <- fs::path(scpca_results, "features", "cell_features_SCPCP000007.tsv")
output_dir <- fs::path(scpca_results, "infercnv")
fs::dir_create(output_dir)
cells_file <- fs::path(output_dir, "scpca_infercnv_cells_SCPCP000007.tsv")
column_summary_file <- fs::path(output_dir, "column_summary.tsv")

stopifnot(
  "Run 08_check-scpca-data.R first" = fs::file_exists(library_summary_file),
  "Run 10_extract-scpca-features.R first" = fs::file_exists(features_file)
)

# Libraries with ScPCA inferCNV results ---------------------------------------------------------
libraries <- data.table::fread(library_summary_file, data.table = FALSE) |>
  dplyr::filter(role == "AML", infercnv_status == "success") |>
  dplyr::mutate(
    file = fs::path(project_dir, sample_id, paste0(library_id, "_processed_rna.h5ad"))
  )
message(nrow(libraries), " AML libraries with ScPCA inferCNV results")
stopifnot(nrow(libraries) == 17, all(fs::file_exists(libraries$file)))

feature_ids <- data.table::fread(features_file, select = "cell_id", data.table = FALSE)$cell_id

# Read an AnnData-encoded data frame (e.g. uns/infercnv_table) with rhdf5 ------------------------
read_h5ad_dataframe <- function(file, group) {
  attrs <- rhdf5::h5readAttributes(file, group)
  stopifnot("Not an AnnData data frame" = identical(attrs[["encoding-type"]], "dataframe"))
  columns <- as.character(attrs[["column-order"]])
  index_name <- as.character(attrs[["_index"]])

  read_column <- function(name) {
    path <- paste0(group, "/", name)
    info <- rhdf5::h5readAttributes(file, path)
    encoding <- info[["encoding-type"]]
    x <- if (identical(encoding, "categorical")) {
      codes <- as.vector(rhdf5::h5read(file, paste0(path, "/codes")))
      categories <- as.vector(rhdf5::h5read(file, paste0(path, "/categories")))
      out <- categories[codes + 1]
      out[codes < 0] <- NA
      out
    } else if (!is.null(encoding) && startsWith(encoding, "nullable")) {
      values <- as.vector(rhdf5::h5read(file, paste0(path, "/values")))
      mask <- as.logical(as.vector(rhdf5::h5read(file, paste0(path, "/mask"))))
      values[mask] <- NA
      values
    } else {
      as.vector(rhdf5::h5read(file, path))
    }
    # HDF5 booleans can come back as factors or strings
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x) && all(x %in% c("TRUE", "FALSE", NA))) x <- as.logical(x)
    x
  }

  df <- as.data.frame(stats::setNames(lapply(columns, read_column), columns), check.names = FALSE)
  if (!"barcodes" %in% columns) {
    df$barcodes <- as.vector(rhdf5::h5read(file, paste0(group, "/", index_name)))
  }
  rhdf5::h5closeAll()
  df
}

# Summarize one column: type, missing values, and its range or most common values
summarize_column <- function(x) {
  if (is.numeric(x)) {
    values <- unique(x[!is.na(x)])
    data.frame(
      type = class(x)[1], n_na = sum(is.na(x)), n_unique = length(values),
      min = suppressWarnings(min(x, na.rm = TRUE)), median = stats::median(x, na.rm = TRUE),
      max = suppressWarnings(max(x, na.rm = TRUE)), frac_nonzero = mean(x != 0, na.rm = TRUE),
      top_values = paste(head(sort(values), 8), collapse = ", ")
    )
  } else {
    tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
    data.frame(
      type = class(x)[1], n_na = sum(is.na(x)), n_unique = length(tab),
      min = NA_real_, median = NA_real_, max = NA_real_, frac_nonzero = NA_real_,
      top_values = paste(head(paste0(names(tab), "=", tab), 8), collapse = ", ")
    )
  }
}

# Extract --------------------------------------------------------------------------------------
cells <- list()
column_summaries <- list()

for (i in seq_len(nrow(libraries))) {
  lib <- libraries[i, ]
  message("[", i, "/", nrow(libraries), "] ", lib$library_id)

  # Cell table: total CNV and reference flag
  adata <- anndataR::read_h5ad(lib$file, as = "HDF5AnnData")
  obs <- adata$obs
  obs_df <- data.frame(
    barcode = adata$obs_names,
    infercnv_total_cnv = obs$infercnv_total_cnv,
    is_infercnv_reference = as.logical(obs$is_infercnv_reference)
  )
  rm(adata, obs)

  # Arm-level table
  arm_table <- read_h5ad_dataframe(lib$file, "uns/infercnv_table")
  keep <- c("subcluster", grep("^(has_(cnv|loss|dupli)|proportion_(loss|dupli))_chr", colnames(arm_table),
                               value = TRUE))
  if (i == 1) {
    message("  uns/infercnv_table: ", nrow(arm_table), " rows, ", ncol(arm_table), " columns; keeping ",
            length(keep))
  }
  arm_df <- arm_table[, c("barcodes", keep)]
  colnames(arm_df)[1] <- "barcode"
  stopifnot("Duplicated barcodes in infercnv_table" = !anyDuplicated(arm_df$barcode))

  n_in_table <- sum(obs_df$barcode %in% arm_df$barcode)
  message("  ", nrow(obs_df), " cells in obs; ", nrow(arm_df), " in infercnv_table; ", n_in_table,
          " in both")
  if (n_in_table < nrow(obs_df)) {
    warning(lib$library_id, ": ", nrow(obs_df) - n_in_table, " cells have no arm-level results")
  }

  df <- obs_df |>
    dplyr::left_join(arm_df, by = "barcode") |>
    dplyr::mutate(cell_id = paste(lib$library_id, barcode, sep = "-"), library_id = lib$library_id,
                  .before = 1)

  n_matched <- sum(df$cell_id %in% feature_ids)
  if (n_matched < nrow(df)) {
    warning(lib$library_id, ": ", nrow(df) - n_matched, " cells not in the feature table")
  }

  summary_cols <- setdiff(colnames(df), c("cell_id", "library_id", "barcode"))
  column_summaries[[lib$library_id]] <- purrr::map(summary_cols, \(col) {
    cbind(data.frame(library_id = lib$library_id, column = col), summarize_column(df[[col]]))
  }) |>
    purrr::list_rbind()

  cells[[lib$library_id]] <- df
}

# Write ----------------------------------------------------------------------------------------
cells <- dplyr::bind_rows(cells)
stopifnot(!anyDuplicated(cells$cell_id))
data.table::fwrite(cells, cells_file, sep = "\t", na = "NA")
message("Wrote ", cells_file, " (", nrow(cells), " cells, ", ncol(cells), " columns)")

column_summaries <- dplyr::bind_rows(column_summaries)
data.table::fwrite(column_summaries, column_summary_file, sep = "\t", na = "NA")
message("Wrote ", column_summary_file)

# What the columns look like: one example of each family, first library
message("\nColumn families (first column of each, first library):")
first_lib <- column_summaries[column_summaries$library_id == libraries$library_id[1], ]
for (prefix in c("infercnv_total_cnv", "is_infercnv_reference", "subcluster", "has_cnv_chr",
                 "has_loss_chr", "has_dupli_chr", "proportion_loss_chr", "proportion_dupli_chr")) {
  row <- first_lib[startsWith(first_lib$column, prefix), ][1, ]
  if (!is.na(row$column)) print(row, row.names = FALSE)
}
