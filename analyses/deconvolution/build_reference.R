#!/usr/bin/env Rscript
# Build a deconvolution reference from labelled single-cell data.
#
# Reads one or more SingleCellExperiment (.rds) or AnnData (.h5ad) files, sums
# raw counts per reference sample x cell state, and saves the reference as a
# SummarizedExperiment (.rds) for deconvolve.R. See SPEC.md.
#
# Example:
#   Rscript build_reference.R \
#     --input lib1_processed_rna.h5ad,lib2_processed_rna.h5ad \
#     --cell-type-column consensus_celltype_annotation \
#     --exclude-labels Unknown \
#     --output references/my_reference.rds

suppressPackageStartupMessages({
  library(optparse)
})

module_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", file_arg)))
})
for (f in c("read_inputs.R", "reference_profiles.R")) {
  source(file.path(module_dir, "R", f))
}

split_list <- function(x) {
  if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",")[[1]])
}

option_list <- list(
  make_option(
    "--input",
    help = "Comma-separated .rds (SingleCellExperiment) or .h5ad files [required]"
  ),
  make_option(
    "--cell-type-column",
    dest = "cell_type_column",
    help = "Cell metadata column holding cell-type labels [required]"
  ),
  make_option(
    "--cell-state-column",
    dest = "cell_state_column",
    help = "Cell metadata column holding cell-state labels [optional]"
  ),
  make_option(
    "--sample-column",
    dest = "sample_column",
    help = paste(
      "Cell metadata column holding reference-sample IDs [optional;",
      "default: a sample_id column, the file's recorded sample ID, or the file name]"
    )
  ),
  make_option(
    "--exclude-labels",
    dest = "exclude_labels",
    default = "",
    help = "Comma-separated cell-type or cell-state labels to drop (e.g. Unknown)"
  ),
  make_option(
    "--unlabelled",
    default = "drop",
    help = "Cells with no cell type: 'drop' or 'unknown' [default: %default]"
  ),
  make_option(
    "--subsample-cap",
    dest = "subsample_cap",
    type = "integer",
    default = 50L,
    help = paste(
      "Max cells per reference sample x cell state in the cell subsample;",
      "0 disables it [default: %default]"
    )
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 2026L,
    help = "Random seed for the subsample [default: %default]"
  ),
  make_option(
    "--name",
    help = "Reference name [default: output file base name]"
  ),
  make_option(
    "--output",
    help = "Output reference .rds [required]"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

for (required in c("input", "cell_type_column", "output")) {
  if (is.null(opts[[required]])) {
    stop("--", gsub("_", "-", required), " is required")
  }
}
if (!opts$unlabelled %in% c("drop", "unknown")) {
  stop("--unlabelled must be 'drop' or 'unknown'")
}
reference_name <- opts$name %||%
  sub("\\.rds$", "", basename(opts$output), ignore.case = TRUE)

reference <- build_reference(
  files = split_list(opts$input),
  cell_type_column = opts$cell_type_column,
  cell_state_column = opts$cell_state_column,
  sample_column = opts$sample_column,
  exclude_labels = split_list(opts$exclude_labels),
  unlabelled = opts$unlabelled,
  subsample_cap = opts$subsample_cap,
  seed = opts$seed,
  reference_name = reference_name
)

dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
saveRDS(reference, opts$output)

groups <- as.data.frame(SummarizedExperiment::colData(reference))
subsample <- S4Vectors::metadata(reference)$cell_subsample
message(
  "Saved reference '",
  reference_name,
  "' to ",
  opts$output,
  ": ",
  nrow(reference),
  " genes, ",
  length(unique(groups$sample_id)),
  " reference samples, ",
  length(unique(groups$cell_type)),
  " cell types, ",
  length(unique(groups$cell_state)),
  " cell states, ",
  sum(groups$n_cells),
  " cells",
  if (!is.null(subsample)) paste0(", subsample of ", ncol(subsample), " cells")
)
