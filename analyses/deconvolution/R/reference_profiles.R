# Building a reference (build_reference.R): summed raw counts per reference
# sample x cell state; and preparing it at run time (deconvolve.R): loading,
# validating, pooling per cell state and applying --min-cells.
#
# The reference is a SummarizedExperiment:
#   assay "counts"  sparse summed counts, genes x (reference sample x cell state)
#   colData         sample_id, cell_type, cell_state, n_cells
#   rowData         gene_id, gene_symbol
#   metadata        deconvolution_reference_version, reference_name, build, cell_subsample

DECONVOLUTION_REFERENCE_VERSION <- 1L

#' Build a reference from one or more labelled single-cell files
#'
#' Files are read one at a time and summed per (reference sample, cell type,
#' cell state) group before the next is read, so peak memory is about one file.
#'
#' @param files Paths to `.rds` (SingleCellExperiment) or `.h5ad` files.
#' @param cell_type_column Cell metadata column with cell-type labels.
#' @param cell_state_column Optional cell metadata column with cell-state labels.
#' @param sample_column Optional cell metadata column with reference-sample IDs.
#' @param exclude_labels Cell-type or cell-state labels whose cells are dropped.
#' @param unlabelled "drop" or "unknown": how to handle cells with no cell type.
#' @param subsample_cap Max cells per reference sample x cell state kept in the
#'   cell subsample; 0 disables the subsample.
#' @param seed Random seed for the subsample.
#' @param reference_name Name stored in the reference.
#' @param small_label_cells Labels with fewer pooled cells are flagged in the
#'   log (nothing is dropped at build time).
#' @return A SummarizedExperiment.
build_reference <- function(
  files,
  cell_type_column,
  cell_state_column = NULL,
  sample_column = NULL,
  exclude_labels = character(0),
  unlabelled = c("drop", "unknown"),
  subsample_cap = 50L,
  seed = 2026L,
  reference_name = "reference",
  small_label_cells = 10L
) {
  unlabelled <- match.arg(unlabelled)
  set.seed(seed)

  sums <- list()
  subsamples <- list()
  genes <- NULL
  n_dropped <- 0L

  for (file in files) {
    message("Reading ", file)
    sc <- read_sc_file(file)

    if (is.null(genes)) {
      genes <- sc$genes
    } else if (!identical(genes$gene_id, sc$genes$gene_id)) {
      stop(
        "Gene IDs in ",
        file,
        " differ from those in ",
        files[1],
        "; all inputs must share the same genes in the same order"
      )
    }

    labels <- resolve_labels(
      sc$cells,
      cell_type_column = cell_type_column,
      cell_state_column = cell_state_column,
      exclude_labels = exclude_labels,
      unlabelled = unlabelled,
      file = file
    )
    labels$sample_id <- resolve_sample_ids(
      sc$cells,
      sample_column = sample_column,
      file_sample_ids = sc$file_sample_ids,
      file_library_id = sc$file_library_id,
      file = file
    )
    keep <- !is.na(labels$cell_type)
    n_dropped <- n_dropped + sum(!keep)
    labels <- labels[keep, , drop = FALSE]
    counts <- sc$counts[, keep, drop = FALSE]
    rm(sc)

    sums[[file]] <- sum_by_group(counts, labels)
    if (subsample_cap > 0) {
      subsamples[[file]] <- draw_subsample(counts, labels, subsample_cap)
    }
    rm(counts)
    gc(verbose = FALSE)
  }

  message("Dropped ", n_dropped, " cells without a usable label")
  summed <- combine_group_sums(sums)
  shared_states <- shared_state_names(
    summed$groups$cell_type,
    summed$groups$cell_state
  )
  summed$groups$cell_state <- disambiguate_states(
    summed$groups$cell_type,
    summed$groups$cell_state,
    shared_states
  )

  cell_subsample <- NULL
  if (subsample_cap > 0) {
    cell_subsample <- combine_subsamples(
      subsamples,
      genes = genes,
      cap = subsample_cap,
      shared_states = shared_states
    )
  }

  flag_small_labels(summed$groups, small_label_cells)

  colnames(summed$counts) <- with(
    summed$groups,
    paste(sample_id, cell_type, cell_state, sep = "|")
  )
  rownames(summed$counts) <- genes$gene_id
  reference <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = summed$counts),
    colData = S4Vectors::DataFrame(
      summed$groups,
      row.names = colnames(summed$counts)
    ),
    rowData = S4Vectors::DataFrame(genes, row.names = genes$gene_id)
  )
  S4Vectors::metadata(reference) <- list(
    deconvolution_reference_version = DECONVOLUTION_REFERENCE_VERSION,
    reference_name = reference_name,
    build = list(
      input_files = normalizePath(files),
      input_md5 = unname(tools::md5sum(files)),
      cell_type_column = cell_type_column,
      cell_state_column = cell_state_column,
      sample_column = sample_column,
      unlabelled = unlabelled,
      exclude_labels = exclude_labels,
      subsample_cap = subsample_cap,
      seed = seed,
      n_cells_dropped = n_dropped,
      package_versions = package_versions(),
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    cell_subsample = cell_subsample
  )
  reference
}

#' Assign each cell a cell type and cell state; NA means the cell is dropped
resolve_labels <- function(
  cells,
  cell_type_column,
  cell_state_column,
  exclude_labels,
  unlabelled,
  file
) {
  missing_columns <- setdiff(
    c(cell_type_column, cell_state_column),
    names(cells)
  )
  if (length(missing_columns) > 0) {
    stop(
      "Column(s) not found in ",
      file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  cell_type <- as.character(cells[[cell_type_column]])
  cell_state <- if (is.null(cell_state_column)) {
    cell_type
  } else {
    as.character(cells[[cell_state_column]])
  }

  unlabelled_cells <- is.na(cell_type) | cell_type == ""
  if (unlabelled == "unknown") {
    cell_type[unlabelled_cells] <- "unknown"
    cell_state[unlabelled_cells] <- "unknown"
  } else {
    cell_type[unlabelled_cells] <- NA
  }
  # a cell type without a state label uses the cell type as its state
  no_state <- !is.na(cell_type) & (is.na(cell_state) | cell_state == "")
  cell_state[no_state] <- cell_type[no_state]

  excluded <- cell_type %in% exclude_labels | cell_state %in% exclude_labels
  cell_type[excluded] <- NA
  cell_state[is.na(cell_type)] <- NA
  data.frame(cell_type = cell_type, cell_state = cell_state)
}

#' Reference-sample ID for each cell
#'
#' Order of precedence: `sample_column`, a `sample_id` cell column, a single
#' sample ID recorded for the file, then the file's base name. A file that
#' records several sample IDs (a multiplexed library) is treated as one
#' reference sample, named by its library ID.
resolve_sample_ids <- function(
  cells,
  sample_column,
  file_sample_ids,
  file_library_id,
  file
) {
  if (!is.null(sample_column)) {
    if (!sample_column %in% names(cells)) {
      stop("Sample column '", sample_column, "' not found in ", file)
    }
    return(as.character(cells[[sample_column]]))
  }
  if ("sample_id" %in% names(cells)) {
    return(as.character(cells$sample_id))
  }
  if (length(file_sample_ids) == 1) {
    return(rep(as.character(file_sample_ids), nrow(cells)))
  }
  file_id <- if (length(file_library_id) == 1) {
    as.character(file_library_id)
  } else {
    sub("\\.(rds|h5ad)$", "", basename(file), ignore.case = TRUE)
  }
  if (length(file_sample_ids) > 1) {
    message(
      file,
      " records several sample IDs (",
      paste(file_sample_ids, collapse = ", "),
      "); treating the library as one reference sample (",
      file_id,
      "). Use --sample-column to split it by sample."
    )
  } else {
    message(
      "No sample ID found in ",
      file,
      "; using ",
      file_id,
      " as its reference sample"
    )
  }
  rep(file_id, nrow(cells))
}

#' Sum counts per (sample_id, cell_type, cell_state) group
#'
#' @return list(counts = genes x groups sparse matrix, groups = data.frame
#'   with sample_id, cell_type, cell_state, n_cells)
sum_by_group <- function(counts, labels) {
  key <- paste(
    labels$sample_id,
    labels$cell_type,
    labels$cell_state,
    sep = "\r"
  )
  group_keys <- unique(key)
  indicator <- Matrix::sparseMatrix(
    i = seq_along(key),
    j = match(key, group_keys),
    x = 1,
    dims = c(length(key), length(group_keys))
  )
  first <- match(group_keys, key)
  list(
    counts = counts %*% indicator,
    groups = data.frame(
      sample_id = labels$sample_id[first],
      cell_type = labels$cell_type[first],
      cell_state = labels$cell_state[first],
      n_cells = as.integer(Matrix::colSums(indicator))
    )
  )
}

#' Combine per-file group sums, summing groups that recur across files
#' (e.g. several libraries from one reference sample)
combine_group_sums <- function(sums) {
  counts <- do.call(cbind, lapply(sums, `[[`, "counts"))
  groups <- do.call(rbind, lapply(sums, `[[`, "groups"))
  key <- paste(
    groups$sample_id,
    groups$cell_type,
    groups$cell_state,
    sep = "\r"
  )
  group_keys <- unique(key)
  if (length(group_keys) == length(key)) {
    rownames(groups) <- NULL
    return(list(counts = counts, groups = groups))
  }
  indicator <- Matrix::sparseMatrix(
    i = seq_along(key),
    j = match(key, group_keys),
    x = 1,
    dims = c(length(key), length(group_keys))
  )
  first <- match(group_keys, key)
  combined <- groups[first, c("sample_id", "cell_type", "cell_state")]
  combined$n_cells <- as.integer(tapply(
    groups$n_cells,
    factor(key, group_keys),
    sum
  ))
  rownames(combined) <- NULL
  list(counts = counts %*% indicator, groups = combined)
}

#' Cell-state names used under more than one cell type
shared_state_names <- function(cell_type, cell_state) {
  pairs <- unique(data.frame(cell_type = cell_type, cell_state = cell_state))
  shared <- unique(pairs$cell_state[duplicated(pairs$cell_state)])
  if (length(shared) > 0) {
    message(
      "Cell-state names shared across cell types will be prefixed with their ",
      "cell type: ",
      paste(shared, collapse = ", ")
    )
  }
  shared
}

#' Make cell-state names unique across cell types
#'
#' States named in `shared` become "cell_type:cell_state"; others are kept.
disambiguate_states <- function(cell_type, cell_state, shared) {
  ifelse(
    cell_state %in% shared,
    paste(cell_type, cell_state, sep = ":"),
    cell_state
  )
}

#' Draw up to `cap` cells per (sample, cell type, cell state) from one file
draw_subsample <- function(counts, labels, cap) {
  key <- paste(
    labels$sample_id,
    labels$cell_type,
    labels$cell_state,
    sep = "\r"
  )
  picked <- unlist(
    lapply(split(seq_along(key), key), function(idx) {
      if (length(idx) > cap) sample(idx, cap) else idx
    }),
    use.names = FALSE
  )
  picked <- sort(picked)
  list(
    counts = counts[, picked, drop = FALSE],
    labels = labels[picked, , drop = FALSE]
  )
}

#' Combine per-file subsamples into one SingleCellExperiment, re-capping groups
#' that span several files, and applying the final cell-state names
combine_subsamples <- function(subsamples, genes, cap, shared_states) {
  counts <- do.call(cbind, lapply(subsamples, `[[`, "counts"))
  labels <- do.call(rbind, lapply(subsamples, `[[`, "labels"))
  key <- paste(
    labels$sample_id,
    labels$cell_type,
    labels$cell_state,
    sep = "\r"
  )
  picked <- sort(unlist(
    lapply(split(seq_along(key), key), function(idx) {
      if (length(idx) > cap) sample(idx, cap) else idx
    }),
    use.names = FALSE
  ))
  counts <- counts[, picked, drop = FALSE]
  labels <- labels[picked, , drop = FALSE]
  # carry over disambiguated state names
  labels$cell_state <- disambiguate_states(
    labels$cell_type,
    labels$cell_state,
    shared_states
  )
  rownames(counts) <- genes$gene_id
  rownames(labels) <- NULL
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(labels),
    rowData = S4Vectors::DataFrame(genes, row.names = genes$gene_id)
  )
}

#' Log cell types and cell states with few pooled cells (nothing is dropped)
flag_small_labels <- function(groups, threshold) {
  by_type <- tapply(groups$n_cells, groups$cell_type, sum)
  # states that are just their cell type (no state column) are covered above
  groups <- groups[groups$cell_state != groups$cell_type, , drop = FALSE]
  by_state <- tapply(
    groups$n_cells,
    paste(groups$cell_type, groups$cell_state, sep = " / "),
    sum
  )
  small_types <- by_type[by_type < threshold]
  small_states <- by_state[by_state < threshold]
  if (length(small_types) > 0) {
    message(
      "Cell types with fewer than ",
      threshold,
      " cells: ",
      paste0(names(small_types), " (", small_types, ")", collapse = ", ")
    )
  }
  if (length(small_states) > 0) {
    message(
      "Cell states with fewer than ",
      threshold,
      " cells: ",
      paste0(names(small_states), " (", small_states, ")", collapse = ", ")
    )
  }
  invisible(list(types = small_types, states = small_states))
}

#' Versions of the packages that shape the reference
package_versions <- function() {
  pkgs <- c(
    "SummarizedExperiment",
    "SingleCellExperiment",
    "anndataR",
    "rhdf5",
    "Matrix"
  )
  versions <- vapply(
    pkgs,
    function(p) as.character(utils::packageVersion(p)),
    character(1)
  )
  c(R = as.character(getRversion()), versions)
}

# Preparing reference profiles at run time (deconvolve.R) ---------------------

REFERENCE_COLDATA_COLUMNS <- c(
  "sample_id",
  "cell_type",
  "cell_state",
  "n_cells"
)

#' Load a reference .rds and check it is one this tool can use
#'
#' @param path Path to a reference written by build_reference.R.
#' @return The reference SummarizedExperiment.
load_reference <- function(path) {
  if (!file.exists(path)) {
    stop("Reference file not found: ", path)
  }
  reference <- readRDS(path)
  validate_reference(reference, path)
  reference
}

#' Stop with a clear message unless `reference` matches the reference format
validate_reference <- function(reference, path) {
  fail <- function(...) {
    stop(
      path,
      " is not a usable deconvolution reference: ",
      ...,
      ". Build references with build_reference.R.",
      call. = FALSE
    )
  }
  if (!methods::is(reference, "SummarizedExperiment")) {
    fail("expected a SummarizedExperiment, found ", class(reference)[1])
  }
  version <- S4Vectors::metadata(reference)$deconvolution_reference_version
  if (is.null(version)) {
    fail("no deconvolution_reference_version in metadata()")
  }
  if (!identical(as.integer(version), DECONVOLUTION_REFERENCE_VERSION)) {
    fail(
      "reference format version ",
      version,
      " is not supported (this tool reads version ",
      DECONVOLUTION_REFERENCE_VERSION,
      ")"
    )
  }
  if (!"counts" %in% SummarizedExperiment::assayNames(reference)) {
    fail("no 'counts' assay")
  }
  missing_col <- setdiff(
    REFERENCE_COLDATA_COLUMNS,
    names(SummarizedExperiment::colData(reference))
  )
  if (length(missing_col) > 0) {
    fail("colData lacks column(s) ", paste(missing_col, collapse = ", "))
  }
  missing_row <- setdiff(
    c("gene_id", "gene_symbol"),
    names(SummarizedExperiment::rowData(reference))
  )
  if (length(missing_row) > 0) {
    fail("rowData lacks column(s) ", paste(missing_row, collapse = ", "))
  }
  if (ncol(reference) == 0) {
    fail("it has no reference columns")
  }
  invisible(reference)
}

#' Whether a reference has cell states distinct from its cell types
has_cell_states <- function(reference) {
  groups <- SummarizedExperiment::colData(reference)
  any(groups$cell_state != groups$cell_type)
}

#' Pool reference columns into one profile per cell state, applying --min-cells
#'
#' Steps, in order: collapse cell states to cell types when `use_cell_states`
#' is FALSE; pool across reference samples; drop cell types, then cell states,
#' with fewer than `min_cells` cells, warning for each level. A cell type whose
#' states are all dropped is dropped too. When states are off (or the
#' reference has none) each cell type has one state named after it.
#'
#' @param reference A reference SummarizedExperiment (see load_reference()).
#' @param use_cell_states Keep cell states (TRUE) or collapse them (FALSE).
#' @param min_cells Minimum pooled cells per cell type and per cell state.
#' @return A list with
#'   `counts`: sparse summed counts, genes x cell states (columns named by
#'     cell state), with the reference's gene IDs as row names;
#'   `labels`: data.frame(cell_type, cell_state, n_cells, n_samples), one row
#'     per column of `counts`;
#'   `genes`: the reference rowData (gene_id, gene_symbol) as a data.frame;
#'   `dropped`: data.frame(level, label, n_cells) of labels removed by
#'     `min_cells`.
prepare_profiles <- function(reference, use_cell_states, min_cells) {
  groups <- as.data.frame(SummarizedExperiment::colData(reference))
  counts <- SummarizedExperiment::assay(reference, "counts")
  if (!use_cell_states) {
    groups$cell_state <- groups$cell_type
  }

  key <- paste(groups$cell_type, groups$cell_state, sep = "\r")
  label_keys <- unique(key)
  indicator <- Matrix::sparseMatrix(
    i = seq_along(key),
    j = match(key, label_keys),
    x = 1,
    dims = c(length(key), length(label_keys))
  )
  first <- match(label_keys, key)
  labels <- data.frame(
    cell_type = groups$cell_type[first],
    cell_state = groups$cell_state[first],
    n_cells = as.integer(tapply(groups$n_cells, factor(key, label_keys), sum)),
    n_samples = as.integer(tapply(
      groups$sample_id,
      factor(key, label_keys),
      function(x) length(unique(x))
    ))
  )
  pooled <- methods::as(counts %*% indicator, "CsparseMatrix")

  kept <- apply_min_cells(labels, min_cells)
  by_label <- order(labels$cell_type[kept$keep], labels$cell_state[kept$keep])
  keep <- which(kept$keep)[by_label]

  pooled <- pooled[, keep, drop = FALSE]
  labels <- labels[keep, , drop = FALSE]
  rownames(labels) <- NULL
  dimnames(pooled) <- list(rownames(reference), labels$cell_state)

  list(
    counts = pooled,
    labels = labels,
    genes = as.data.frame(SummarizedExperiment::rowData(reference))[, c(
      "gene_id",
      "gene_symbol"
    )],
    dropped = kept$dropped
  )
}

#' Decide which pooled labels pass --min-cells, warning about the rest
#'
#' @param labels data.frame(cell_type, cell_state, n_cells), one row per state.
#' @return list(keep = logical per row, dropped = data.frame(level, label,
#'   n_cells)).
apply_min_cells <- function(labels, min_cells) {
  type_cells <- tapply(labels$n_cells, labels$cell_type, sum)
  small_types <- type_cells[type_cells < min_cells]
  keep <- !labels$cell_type %in% names(small_types)

  # states are checked separately only where they differ from their type
  has_state <- labels$cell_state != labels$cell_type
  small_state <- keep & has_state & labels$n_cells < min_cells
  keep <- keep & !small_state

  # a type that passed but lost every state is gone too
  emptied <- setdiff(
    unique(labels$cell_type[!labels$cell_type %in% names(small_types)]),
    labels$cell_type[keep]
  )

  dropped <- rbind(
    data.frame(
      level = rep("cell_type", length(small_types)),
      label = names(small_types),
      n_cells = as.integer(small_types)
    ),
    data.frame(
      level = rep("cell_state", sum(small_state)),
      label = labels$cell_state[small_state],
      n_cells = labels$n_cells[small_state]
    ),
    data.frame(
      level = rep("cell_type", length(emptied)),
      label = emptied,
      n_cells = as.integer(type_cells[emptied])
    )
  )

  if (length(small_types) > 0) {
    warning(
      "Dropping cell type(s) with fewer than ",
      min_cells,
      " cells (--min-cells): ",
      paste0(names(small_types), " (", small_types, ")", collapse = ", "),
      call. = FALSE
    )
  }
  if (any(small_state)) {
    warning(
      "Dropping cell state(s) with fewer than ",
      min_cells,
      " cells (--min-cells): ",
      paste0(
        labels$cell_state[small_state],
        " (",
        labels$n_cells[small_state],
        ")",
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  if (length(emptied) > 0) {
    warning(
      "Cell type(s) left with no cell states after --min-cells, so dropped: ",
      paste(emptied, collapse = ", "),
      ". Consider --no-cell-states.",
      call. = FALSE
    )
  }
  if (!any(keep)) {
    stop("No cell types left after applying --min-cells ", min_cells)
  }
  list(keep = keep, dropped = dropped)
}
