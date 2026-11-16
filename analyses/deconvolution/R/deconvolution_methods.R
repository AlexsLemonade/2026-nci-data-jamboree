# Deconvolution methods and the method registry (SPEC.md §6 step 6, §7).
#
# Every method has the signature run_<method>(reference, bulk, options):
#   reference  list(counts, labels, markers):
#                counts   summed counts, genes x labels (sparse or dense),
#                         columns named by cell state (= cell type when cell
#                         states are off), after gene harmonisation/exclusion
#                labels   data.frame(cell_type, cell_state, ...), one row per
#                         column of `counts` (see prepare_profiles())
#                markers  marker genes, for methods with needs_markers = TRUE
#   bulk       numeric matrix, genes x bulk samples, same rows as counts
#   options    the method's own options, named as in the registry
# and returns list(
#   proportions = data.frame(method, level, label, sample_id, quantity,
#                            proportion),
#   info        = list of method details for run info (may be empty),
#   object      = the method's raw result, for --save-method-objects
# ). deconvolve.R adds the `reference` column to `proportions`.

#' The method registry
#'
#' For each method: the function that runs it, its options and defaults, the
#' quantity it estimates, whether it needs marker genes or the cell
#' subsample, and the levels it reports. Methods reporting only `cell_type`
#' are given references with cell states collapsed.
DECONVOLUTION_METHODS <- list(
  instaprism = list(
    run = function(...) run_instaprism(...),
    options = list(
      instaprism_filter = TRUE,
      malignant_cell_type = NULL,
      update_mode = "all",
      n_cores = 1L
    ),
    quantity = "mrna_share",
    needs_markers = FALSE,
    needs_cell_subsample = FALSE,
    levels = c("cell_type", "cell_state")
  ),
  nnls = list(
    run = function(...) run_nnls(...),
    options = list(),
    quantity = "mrna_share",
    needs_markers = TRUE,
    needs_cell_subsample = FALSE,
    levels = "cell_type"
  )
)

#' Look up a method in the registry, failing clearly for unknown names
get_method <- function(name) {
  if (!name %in% names(DECONVOLUTION_METHODS)) {
    stop(
      "Unknown method '",
      name,
      "'. Available: ",
      paste(names(DECONVOLUTION_METHODS), collapse = ", "),
      call. = FALSE
    )
  }
  DECONVOLUTION_METHODS[[name]]
}

#' Whether a method reports cell-state proportions
method_uses_cell_states <- function(name) {
  "cell_state" %in% get_method(name)$levels
}

#' Fill a method's options from its registry defaults
#'
#' @param name Method name.
#' @param options Named list of option values; NULL entries and options the
#'   method does not declare are ignored.
#' @return The method's options, with defaults for anything not given.
method_options <- function(name, options) {
  defaults <- get_method(name)$options
  given <- options[intersect(names(options), names(defaults))]
  given <- given[!vapply(given, is.null, logical(1))]
  utils::modifyList(defaults, given)
}

#' Run a method by name with its options
run_method <- function(name, reference, bulk, options) {
  method <- get_method(name)
  if (method$needs_markers && length(reference$markers) == 0) {
    stop("Method '", name, "' needs marker genes", call. = FALSE)
  }
  method$run(reference, bulk, method_options(name, options))
}

#' Convert a labels x samples proportion matrix to the long table
theta_to_long <- function(theta, method, level, quantity) {
  theta <- as.matrix(theta)
  data.frame(
    method = method,
    level = level,
    label = rep(rownames(theta), times = ncol(theta)),
    sample_id = rep(colnames(theta), each = nrow(theta)),
    quantity = quantity,
    proportion = as.vector(theta)
  )
}

# InstaPrism -------------------------------------------------------------------

#' Deconvolve with InstaPrism, with the optional reference update
#'
#' The reference is built with refPrepare() from summed counts per cell
#' state. The initial fit is reported as method `instaprism` (cell_type, plus
#' cell_state when the reference has cell states). With a malignant cell
#' type, InstaPrism_update() re-estimates the malignant profile per bulk
#' sample (and, in `all` mode, the normal profiles pooled across the run);
#' it is reported as `instaprism_updated` (see SPEC.md §6 step 6).
#'
#' Options: `instaprism_filter`, `malignant_cell_type`, `update_mode`
#' (`all` or `malignant`), `n_cores`.
run_instaprism <- function(reference, bulk, options) {
  suppressPackageStartupMessages(library(InstaPrism))
  quantity <- DECONVOLUTION_METHODS$instaprism$quantity
  labels <- reference$labels
  use_states <- any(labels$cell_state != labels$cell_type)
  malignant <- options$malignant_cell_type
  if (!is.null(malignant) && !malignant %in% labels$cell_type) {
    stop(
      "--malignant-cell-type '",
      malignant,
      "' is not a cell type in the reference. Cell types: ",
      paste(unique(labels$cell_type), collapse = ", "),
      call. = FALSE
    )
  }
  if (!options$update_mode %in% c("all", "malignant")) {
    stop("--update-mode must be 'all' or 'malignant'", call. = FALSE)
  }
  if (options$instaprism_filter && ncol(bulk) < 10) {
    warning(
      "InstaPrism's filter is on with only ",
      ncol(bulk),
      " bulk samples: a gene above 1% of reads in a single sample is ",
      "dropped for every sample. Consider --no-instaprism-filter.",
      call. = FALSE
    )
  }

  ref_counts <- as.matrix(reference$counts)
  ref_phi <- InstaPrism::refPrepare(
    sc_Expr = ref_counts,
    cell.type.labels = labels$cell_type,
    cell.state.labels = labels$cell_state
  )
  bulk <- as.matrix(bulk)
  # InstaPrism forces text progress bars on stdout; keep them out of the log
  fit <- quietly(InstaPrism::InstaPrism(
    bulk_Expr = bulk,
    refPhi_cs = ref_phi,
    filter = options$instaprism_filter,
    n.core = options$n_cores,
    verbose = FALSE
  ))

  proportions <- list(theta_to_long(
    fit@Post.ini.ct@theta,
    "instaprism",
    "cell_type",
    quantity
  ))
  if (use_states) {
    proportions <- c(
      proportions,
      list(theta_to_long(
        fit@Post.ini.cs@theta,
        "instaprism",
        "cell_state",
        quantity
      ))
    )
  }

  # genes refPrepare() drops (zero in the reference) and genes the filter drops
  reference_genes <- rownames(ref_counts)[rowSums(ref_counts) > 0]
  used_genes <- rownames(fit@initial.scaler)
  info <- list(
    n_genes_used = length(used_genes),
    genes_zero_in_reference = setdiff(rownames(ref_counts), reference_genes),
    genes_removed_by_filter = setdiff(
      intersect(reference_genes, rownames(bulk)),
      used_genes
    ),
    filter = options$instaprism_filter,
    malignant_cell_type = malignant,
    update_mode = if (!is.null(malignant)) options$update_mode
  )
  object <- list(initial = fit)

  if (!is.null(malignant)) {
    # the update matches bulk samples by column position: pass the same matrix
    updated <- quietly(InstaPrism::InstaPrism_update(
      InstaPrism_obj = fit,
      bulk_Expr = bulk,
      key = malignant,
      n.core = options$n_cores,
      cell.types.to.update = if (options$update_mode == "all") "all" else NULL
    ))
    object$updated <- updated
    proportions <- c(
      proportions,
      list(theta_to_long(
        InstaPrism::merge_updated_theta(updated),
        "instaprism_updated",
        "cell_type",
        quantity
      ))
    )
    if (options$update_mode == "malignant" && use_states) {
      # normal cell types keep their cell states; look rows up by name
      normal_states <- labels$cell_state[labels$cell_type != malignant]
      state_rows <- intersect(normal_states, rownames(updated@theta))
      proportions <- c(
        proportions,
        list(theta_to_long(
          updated@theta[state_rows, , drop = FALSE],
          "instaprism_updated",
          "cell_state",
          quantity
        ))
      )
    }
  }

  list(
    proportions = do.call(rbind, proportions),
    info = info,
    object = object
  )
}

#' Evaluate an expression, discarding what it prints to stdout
quietly <- function(expr) {
  utils::capture.output(result <- expr)
  result
}

# NNLS -------------------------------------------------------------------------

#' Deconvolve with non-negative least squares on marker genes
#'
#' Profiles are the CPM of the summed counts per cell type, and each bulk
#' sample is CPM-normalised over the same genes; both are then restricted to
#' the markers. Coefficients from nnls::nnls() are renormalised to sum to 1.
#' Fits at cell-type level: `reference` should have cell states collapsed.
run_nnls <- function(reference, bulk, options) {
  quantity <- DECONVOLUTION_METHODS$nnls$quantity
  markers <- reference$markers
  profiles <- to_cpm(as.matrix(reference$counts))[markers, , drop = FALSE]
  colnames(profiles) <- reference$labels$cell_type
  bulk_cpm <- to_cpm(as.matrix(bulk))[markers, , drop = FALSE]

  theta <- apply(bulk_cpm, 2, function(y) {
    coef <- nnls::nnls(profiles, y)$x
    if (sum(coef) > 0) coef / sum(coef) else rep(NA_real_, length(coef))
  })
  theta <- matrix(
    theta,
    nrow = ncol(profiles),
    dimnames = list(colnames(profiles), colnames(bulk))
  )
  failed <- colnames(theta)[is.na(theta[1, ])]
  if (length(failed) > 0) {
    warning(
      "NNLS gave all-zero coefficients for bulk sample(s): ",
      paste(utils::head(failed, 5), collapse = ", "),
      "; their proportions are NA",
      call. = FALSE
    )
  }

  list(
    proportions = theta_to_long(theta, "nnls", "cell_type", quantity),
    info = list(n_markers_used = length(markers)),
    object = theta
  )
}

#' Counts per million, per column
to_cpm <- function(counts) {
  sweep(counts, 2, colSums(counts), "/") * 1e6
}
