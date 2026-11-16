# Gene selection for deconvolve.R: matching bulk and reference genes, gene
# category exclusion, and marker selection (SPEC.md §6 steps 4-5, §7, §8).
#
# The harmonised genes object passed between these functions is a list:
#   reference  sparse summed counts, genes x labels (from prepare_profiles())
#   bulk       numeric matrix, genes x bulk samples, same rows as `reference`
#   genes      data.frame(gene, gene_symbol): `gene` is the shared key
#              (Ensembl ID or symbol, per gene_id_type), one row per row above
#   gene_id_type  "ensembl" or "symbol"

#' Match bulk genes to reference genes and keep the shared genes
#'
#' Bulk IDs are matched against the reference's `gene_id` (Ensembl) or
#' `gene_symbol`, by the bulk's detected ID type. When matching by symbol,
#' reference genes without a symbol are dropped and reference genes sharing a
#' symbol are summed, with a warning. Genes are kept in bulk order.
#'
#' @param profiles Output of prepare_profiles().
#' @param bulk Output of read_bulk_matrix().
#' @param min_gene_overlap Fail if fewer than this fraction of bulk genes
#'   match the reference.
#' @return A harmonised genes object (see top of file) plus `funnel`, a list
#'   of gene counts: `bulk`, `reference`, `shared`.
harmonise_genes <- function(profiles, bulk, min_gene_overlap) {
  reference <- profiles$counts
  gene_id_type <- bulk$gene_id_type
  if (gene_id_type == "ensembl") {
    symbols <- stats::setNames(
      profiles$genes$gene_symbol,
      profiles$genes$gene_id
    )
  } else {
    reference <- sum_reference_by_symbol(reference, profiles$genes$gene_symbol)
    symbols <- stats::setNames(rownames(reference), rownames(reference))
  }

  bulk_counts <- bulk$counts
  shared <- rownames(bulk_counts)[
    rownames(bulk_counts) %in% rownames(reference)
  ]
  overlap <- length(shared) / nrow(bulk_counts)
  message(
    length(shared),
    " of ",
    nrow(bulk_counts),
    " bulk genes (",
    sprintf("%.1f%%", 100 * overlap),
    ") match the reference by ",
    if (gene_id_type == "ensembl") "Ensembl ID" else "symbol",
    "; ",
    nrow(bulk_counts) - length(shared),
    " bulk and ",
    nrow(reference) - length(shared),
    " reference genes are not shared"
  )
  if (overlap < min_gene_overlap) {
    stop(
      "Only ",
      sprintf("%.1f%%", 100 * overlap),
      " of bulk genes match the reference (--min-gene-overlap ",
      min_gene_overlap,
      "). Check that the bulk and reference use the same gene IDs and species.",
      call. = FALSE
    )
  }

  bulk_counts <- bulk_counts[shared, , drop = FALSE]
  check_bulk_depth(bulk_counts)
  list(
    reference = reference[shared, , drop = FALSE],
    bulk = bulk_counts,
    genes = data.frame(gene = shared, gene_symbol = unname(symbols[shared])),
    gene_id_type = gene_id_type,
    funnel = list(
      bulk = nrow(bulk$counts),
      reference = nrow(reference),
      shared = length(shared)
    )
  )
}

#' Re-key reference counts by symbol: drop NA symbols, sum duplicates
sum_reference_by_symbol <- function(counts, gene_symbol) {
  has_symbol <- !is.na(gene_symbol) & gene_symbol != ""
  counts <- counts[has_symbol, , drop = FALSE]
  gene_symbol <- gene_symbol[has_symbol]
  message(
    "Matching by symbol: dropped ",
    sum(!has_symbol),
    " reference genes without a symbol"
  )
  duplicated_symbols <- unique(gene_symbol[duplicated(gene_symbol)])
  symbols <- unique(gene_symbol)
  if (length(duplicated_symbols) > 0) {
    warning(
      length(duplicated_symbols),
      " gene symbol(s) are shared by several reference genes; summing their ",
      "counts. Examples: ",
      paste(utils::head(duplicated_symbols, 5), collapse = ", "),
      call. = FALSE
    )
    indicator <- Matrix::sparseMatrix(
      i = match(gene_symbol, symbols),
      j = seq_along(gene_symbol),
      x = 1,
      dims = c(length(symbols), length(gene_symbol))
    )
    counts <- methods::as(indicator %*% counts, "CsparseMatrix")
  }
  rownames(counts) <- symbols
  counts
}

#' Stop if any bulk sample has no counts on the shared genes
check_bulk_depth <- function(bulk_counts) {
  empty <- colnames(bulk_counts)[colSums(bulk_counts) == 0]
  if (length(empty) > 0) {
    stop(
      "Bulk sample(s) with no counts on genes shared with the reference: ",
      paste(utils::head(empty, 5), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(bulk_counts)
}

#' Read a gene category table (gene_id, gene_symbol, category)
read_gene_categories <- function(path) {
  if (!file.exists(path)) {
    stop("Gene category table not found: ", path)
  }
  table <- data.table::fread(
    path,
    sep = "\t",
    data.table = FALSE,
    colClasses = "character"
  )
  missing <- setdiff(c("gene_id", "gene_symbol", "category"), names(table))
  if (length(missing) > 0) {
    stop(
      "Gene category table ",
      path,
      " lacks column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  table$gene_id <- strip_version(table$gene_id)
  table
}

#' Remove genes in the chosen categories
#'
#' Genes are matched to the category table on Ensembl ID when the bulk is
#' keyed by Ensembl ID, and on symbol when it is keyed by symbol.
#'
#' @param harmonised Output of harmonise_genes().
#' @param categories Gene category table (see read_gene_categories()).
#' @param exclude Category names to exclude; `character(0)` excludes nothing.
#' @return `harmonised` with the genes removed, plus `excluded_genes`, a
#'   data.frame(gene, category) with one row per removed gene x matching
#'   category, and `funnel$after_exclusion`.
exclude_gene_categories <- function(harmonised, categories, exclude) {
  unknown <- setdiff(exclude, categories$category)
  if (length(unknown) > 0) {
    stop(
      "Unknown gene categor(ies): ",
      paste(unknown, collapse = ", "),
      ". The gene category table has: ",
      paste(sort(unique(categories$category)), collapse = ", "),
      call. = FALSE
    )
  }
  key <- if (harmonised$gene_id_type == "ensembl") "gene_id" else "gene_symbol"
  chosen <- categories[
    categories$category %in%
      exclude &
      categories[[key]] %in% harmonised$genes$gene,
    ,
    drop = FALSE
  ]
  excluded_genes <- unique(data.frame(
    gene = chosen[[key]],
    category = chosen$category
  ))
  rownames(excluded_genes) <- NULL

  keep <- !harmonised$genes$gene %in% excluded_genes$gene
  if (length(exclude) > 0) {
    by_category <- table(excluded_genes$category)
    message(
      "Excluded ",
      sum(!keep),
      " genes by category: ",
      paste0(names(by_category), " (", by_category, ")", collapse = ", ")
    )
  }
  if (!any(keep)) {
    stop("No genes left after gene category exclusion", call. = FALSE)
  }

  harmonised$reference <- harmonised$reference[keep, , drop = FALSE]
  harmonised$bulk <- harmonised$bulk[keep, , drop = FALSE]
  harmonised$genes <- harmonised$genes[keep, , drop = FALSE]
  rownames(harmonised$genes) <- NULL
  harmonised$excluded_genes <- excluded_genes
  harmonised$funnel$after_exclusion <- sum(keep)
  check_bulk_depth(harmonised$bulk)
  harmonised
}

#' Select marker genes from summed reference profiles
#'
#' For each label (column), genes are ranked by log2(CPM + 1) in that label
#' minus the maximum log2(CPM + 1) over the other labels. Of genes with CPM
#' above `min_cpm` in the label, the top `n_markers` are taken, keeping those
#' whose difference exceeds `min_log2_diff`. Returns the union over labels.
#' CPM is computed over the genes passed in, so call this after gene
#' category exclusion.
#'
#' @param counts Summed counts, genes x labels (at the level being fitted).
#' @return Marker gene names (row names of `counts`), in row order.
select_markers <- function(
  counts,
  n_markers,
  min_cpm = 20,
  min_log2_diff = 0.5
) {
  if (ncol(counts) < 2) {
    stop("Marker selection needs at least two labels", call. = FALSE)
  }
  counts <- as.matrix(counts)
  cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
  log_cpm <- log2(cpm + 1)

  markers <- character(0)
  for (j in seq_len(ncol(counts))) {
    diff <- log_cpm[, j] - matrixStats::rowMaxs(log_cpm[, -j, drop = FALSE])
    diff[cpm[, j] <= min_cpm] <- -Inf
    top <- utils::head(order(diff, decreasing = TRUE), n_markers)
    top <- top[diff[top] > min_log2_diff]
    if (length(top) == 0) {
      warning(
        "No marker genes found for ",
        colnames(counts)[j],
        call. = FALSE
      )
    }
    markers <- union(markers, rownames(counts)[top])
  }
  rownames(counts)[rownames(counts) %in% markers]
}
