# Reading single-cell inputs (SCE .rds and AnnData .h5ad).
#
# Every reader returns the same shape:
#   list(
#     counts  = dgCMatrix of raw integer counts, genes x cells,
#     genes   = data.frame(gene_id, gene_symbol), one row per counts row,
#     cells   = data.frame of per-cell metadata, one row per counts column,
#     file_sample_ids = sample ID(s) recorded for the whole file (may be NULL)
#     file_library_id = library ID recorded for the whole file (may be NULL)
#   )

#' Read raw counts and metadata from one single-cell file
#'
#' @param path Path to a SingleCellExperiment `.rds` or an AnnData `.h5ad`.
#' @return A list with `counts`, `genes`, `cells`, `file_sample_ids` and
#'   `file_library_id`.
read_sc_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  sc <- switch(
    ext,
    rds = read_sce_rds(path),
    h5ad = read_h5ad_counts(path),
    stop("Unsupported single-cell input (expected .rds or .h5ad): ", path)
  )
  sc$counts <- round_counts(sc$counts, path)
  sc$genes$gene_id <- strip_version(sc$genes$gene_id)
  rownames(sc$counts) <- sc$genes$gene_id
  sc
}

#' Read a SingleCellExperiment saved as .rds
read_sce_rds <- function(path) {
  sce <- readRDS(path)
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("Not a SingleCellExperiment: ", path)
  }
  if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
    stop("No 'counts' assay in ", path)
  }
  counts <- methods::as(
    SummarizedExperiment::assay(sce, "counts"),
    "CsparseMatrix"
  )
  row_data <- SummarizedExperiment::rowData(sce)
  gene_id <- if ("gene_ids" %in% names(row_data)) {
    row_data$gene_ids
  } else {
    rownames(sce)
  }
  gene_symbol <- if ("gene_symbol" %in% names(row_data)) {
    as.character(row_data$gene_symbol)
  } else {
    rep(NA_character_, nrow(sce))
  }
  # ScPCA stores sample/library IDs in metadata(), not colData()
  file_meta <- S4Vectors::metadata(sce)
  list(
    counts = counts,
    genes = data.frame(
      gene_id = as.character(gene_id),
      gene_symbol = gene_symbol
    ),
    cells = as.data.frame(SummarizedExperiment::colData(sce)),
    file_sample_ids = file_meta$sample_id,
    file_library_id = file_meta$library_id
  )
}

#' Read raw counts from an AnnData .h5ad
#'
#' Cell metadata comes from anndataR (read lazily, so X is not loaded). Raw
#' counts and gene metadata come from `raw/X` and `raw/var` via rhdf5, because
#' anndataR does not read `.raw`. Files without `raw/X` are read entirely with
#' anndataR, taking counts from `X`.
read_h5ad_counts <- function(path) {
  contents <- rhdf5::h5ls(path, recursive = 2)
  rhdf5::h5closeAll()
  has_raw_x <- any(contents$group == "/raw" & contents$name == "X")

  adata <- anndataR::read_h5ad(path, as = "HDF5AnnData")
  cells <- as.data.frame(adata$obs)
  if (has_raw_x) {
    # release anndataR's file handle before reading with rhdf5
    adata$close()
    counts <- read_h5_sparse(path, "raw/X")
    gene_id <- read_h5_index(path, "raw/var")
    gene_symbol <- read_h5_column(path, "raw/var", "gene_symbol")
    rhdf5::h5closeAll()
  } else {
    counts <- methods::as(Matrix::t(adata$X), "CsparseMatrix")
    var <- as.data.frame(adata$var)
    adata$close()
    gene_id <- if ("gene_ids" %in% names(var)) var$gene_ids else rownames(var)
    gene_symbol <- if ("gene_symbol" %in% names(var)) {
      as.character(var$gene_symbol)
    } else {
      rep(NA_character_, nrow(var))
    }
  }

  if (ncol(counts) != nrow(cells) || nrow(counts) != length(gene_id)) {
    stop(
      "Dimension mismatch in ",
      path,
      ": counts are ",
      nrow(counts),
      " x ",
      ncol(counts),
      " but there are ",
      length(gene_id),
      " genes and ",
      nrow(cells),
      " cells"
    )
  }
  colnames(counts) <- rownames(cells)
  list(
    counts = counts,
    genes = data.frame(
      gene_id = as.character(gene_id),
      gene_symbol = gene_symbol
    ),
    cells = cells,
    file_sample_ids = NULL,
    file_library_id = NULL
  )
}

#' Read a sparse matrix stored in AnnData's CSR/CSC encoding, as genes x cells
read_h5_sparse <- function(path, group) {
  attrs <- rhdf5::h5readAttributes(path, group)
  shape <- as.integer(attrs$shape) # (n_cells, n_genes)
  data <- as.numeric(rhdf5::h5read(path, paste0(group, "/data")))
  indices <- as.integer(rhdf5::h5read(path, paste0(group, "/indices")))
  indptr <- as.numeric(rhdf5::h5read(path, paste0(group, "/indptr")))
  encoding <- attrs[["encoding-type"]]
  if (identical(encoding, "csc_matrix")) {
    # columns of cells x genes are genes; indptr runs over genes
    Matrix::sparseMatrix(
      i = rep(seq_len(shape[2]), diff(indptr)),
      j = indices + 1L,
      x = data,
      dims = c(shape[2], shape[1])
    )
  } else if (identical(encoding, "csr_matrix")) {
    # rows of cells x genes are cells; indptr runs over cells
    Matrix::sparseMatrix(
      i = indices + 1L,
      j = rep(seq_len(shape[1]), diff(indptr)),
      x = data,
      dims = c(shape[2], shape[1])
    )
  } else {
    stop("Unsupported encoding '", encoding, "' for ", group, " in ", path)
  }
}

#' Read the index (row names) of an AnnData dataframe group
read_h5_index <- function(path, group) {
  index_name <- rhdf5::h5readAttributes(path, group)[["_index"]]
  as.character(rhdf5::h5read(path, paste0(group, "/", index_name)))
}

#' Read one column of an AnnData dataframe group, decoding categoricals
#'
#' Returns NA for a missing column. Categorical codes of -1 mean NA.
read_h5_column <- function(path, group, column) {
  column_path <- paste0(group, "/", column)
  contents <- rhdf5::h5ls(path, recursive = TRUE)
  full_paths <- sub("^/", "", file.path(contents$group, contents$name))
  if (!column_path %in% full_paths) {
    n <- length(read_h5_index(path, group))
    return(rep(NA_character_, n))
  }
  attrs <- rhdf5::h5readAttributes(path, column_path)
  if (identical(attrs[["encoding-type"]], "categorical")) {
    categories <- as.character(rhdf5::h5read(
      path,
      paste0(column_path, "/categories")
    ))
    codes <- as.integer(rhdf5::h5read(path, paste0(column_path, "/codes")))
    out <- rep(NA_character_, length(codes))
    out[codes >= 0] <- categories[codes[codes >= 0] + 1L]
    out
  } else {
    as.character(rhdf5::h5read(path, column_path))
  }
}

#' Round non-integer counts, with a warning
#'
#' Raw counts are expected. Non-integer values (e.g. estimated counts) are
#' rounded; the warning reports how many values were affected, since a large
#' share can mean a normalised matrix was read instead of raw counts.
#'
#' @return The counts, rounded if needed.
round_counts <- function(counts, path) {
  values <- if (methods::is(counts, "sparseMatrix")) counts@x else counts
  non_integer <- values != round(values)
  if (any(non_integer)) {
    warning(
      "Counts in ",
      path,
      " are not all integers (",
      sprintf("%.1f%%", 100 * mean(non_integer)),
      " of non-zero values); rounding. Raw counts are expected: check that ",
      "a normalised matrix was not read instead.",
      call. = FALSE
    )
    if (methods::is(counts, "sparseMatrix")) {
      counts@x <- round(counts@x)
      counts <- Matrix::drop0(counts)
    } else {
      counts <- round(counts)
    }
  }
  counts
}

#' Strip Ensembl version suffixes (ENSG00000141510.17 -> ENSG00000141510)
strip_version <- function(ids) {
  sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", ids)
}
