#!/usr/bin/env Rscript
# Estimate cell type and cell state proportions in bulk samples from a
# reference built by build_reference.R.
#
# Writes <prefix>_proportions.tsv and <prefix>_run_info.json to --output-dir,
# plus the HTML run report unless --no-report is given.
#
# Example:
#   Rscript deconvolve.R \
#     --reference references/my_reference.rds \
#     --bulk bulk_counts.tsv \
#     --malignant-cell-type "Malignant" \
#     --output-dir results/

suppressPackageStartupMessages({
  library(optparse)
  # load quietly up front so package start-up warnings don't count as run
  # warnings
  suppressWarnings(requireNamespace("SummarizedExperiment"))
})

module_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", file_arg)))
})
for (f in list.files(file.path(module_dir, "R"), "\\.R$", full.names = TRUE)) {
  source(f)
}

split_list <- function(x) {
  if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",")[[1]])
}

option_list <- list(
  make_option(
    "--reference",
    help = "Reference .rds from build_reference.R [required]"
  ),
  make_option(
    "--bulk",
    help = paste(
      "Genes x bulk samples count matrix: .tsv or .csv (first column = gene",
      "ID) or .rds (matrix/data.frame with gene row names) [required]"
    )
  ),
  make_option(
    "--output-dir",
    dest = "output_dir",
    help = "Directory for output files [required]"
  ),
  make_option(
    "--output-prefix",
    dest = "output_prefix",
    help = "Prefix for output files [default: <bulk>_<reference>_<method>]"
  ),
  make_option(
    "--method",
    default = "instaprism",
    help = paste0(
      "Deconvolution method: ",
      paste(names(DECONVOLUTION_METHODS), collapse = " or "),
      " [default: %default]"
    )
  ),
  make_option(
    "--use-cell-states",
    dest = "use_cell_states",
    action = "store_true",
    default = NA,
    help = paste(
      "Fit cell states [default: on if the reference has cell states and",
      "the method reports them]"
    )
  ),
  make_option(
    "--no-cell-states",
    dest = "use_cell_states",
    action = "store_false",
    help = "Collapse cell states to cell types before fitting"
  ),
  make_option(
    "--min-cells",
    dest = "min_cells",
    type = "integer",
    default = 10L,
    help = "Drop cell types/states with fewer cells [default: %default]"
  ),
  make_option(
    "--malignant-cell-type",
    dest = "malignant_cell_type",
    help = paste(
      "Malignant cell type (a label in the reference). Turns on the",
      "InstaPrism reference update (modes 'all' and 'malignant')"
    )
  ),
  make_option(
    "--update-mode",
    dest = "update_mode",
    help = paste(
      "InstaPrism reference update. With --malignant-cell-type: 'malignant'",
      "updates only the malignant profile, per bulk sample; 'all' also",
      "re-estimates normal profiles, pooled across the run [default].",
      "Without it: 'pooled' re-estimates every profile, pooled across the run"
    )
  ),
  make_option(
    "--instaprism-filter",
    dest = "instaprism_filter",
    action = "store_true",
    default = TRUE,
    help = "Use InstaPrism's outlier-gene filter [default: on]"
  ),
  make_option(
    "--no-instaprism-filter",
    dest = "instaprism_filter",
    action = "store_false",
    help = "Turn off InstaPrism's outlier-gene filter"
  ),
  make_option(
    "--exclude-gene-categories",
    dest = "exclude_gene_categories",
    default = "mitochondrial,ribosomal,sex",
    help = paste(
      "Comma-separated gene categories to exclude; 'none' disables",
      "[default: %default]"
    )
  ),
  make_option(
    "--gene-categories",
    dest = "gene_categories",
    default = file.path(module_dir, "reference-files", "gene_categories.tsv"),
    help = "Gene category table [default: the module's gene_categories.tsv]"
  ),
  make_option(
    "--min-gene-overlap",
    dest = "min_gene_overlap",
    type = "double",
    default = 0.5,
    help = paste(
      "Fail if fewer than this fraction of bulk genes match the reference",
      "[default: %default]"
    )
  ),
  make_option(
    "--n-markers",
    dest = "n_markers",
    type = "integer",
    default = 60L,
    help = "Markers per label, for methods that use markers [default: %default]"
  ),
  make_option(
    "--n-cores",
    dest = "n_cores",
    type = "integer",
    default = 1L,
    help = "Cores for InstaPrism; does not change results [default: %default]"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 2026L,
    help = "Random seed [default: %default]"
  ),
  make_option(
    "--save-method-objects",
    dest = "save_method_objects",
    action = "store_true",
    default = FALSE,
    help = "Also save the method's raw result object as .rds"
  ),
  make_option(
    "--no-report",
    dest = "no_report",
    action = "store_true",
    default = FALSE,
    help = "Skip the HTML run report"
  ),
  make_option(
    "--save-figures",
    dest = "save_figures",
    action = "store_true",
    default = FALSE,
    help = "Also save each report figure as a separate file"
  ),
  make_option(
    "--figure-format",
    dest = "figure_format",
    default = "png",
    help = "Figure format for --save-figures: png or pdf [default: %default]"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))
options(warn = 1) # print warnings as they happen

# Option checks ----------------------------------------------------------------

for (required in c("reference", "bulk", "output_dir")) {
  if (is.null(opts[[required]])) {
    stop("--", gsub("_", "-", required), " is required", call. = FALSE)
  }
}
method <- get_method(opts$method)
if (
  !is.null(opts$update_mode) &&
    !opts$update_mode %in% c("all", "malignant", "pooled")
) {
  stop("--update-mode must be 'all', 'malignant' or 'pooled'", call. = FALSE)
}
if (!opts$figure_format %in% c("png", "pdf")) {
  stop("--figure-format must be 'png' or 'pdf'", call. = FALSE)
}
if (opts$min_gene_overlap < 0 || opts$min_gene_overlap > 1) {
  stop("--min-gene-overlap must be between 0 and 1", call. = FALSE)
}
exclude_categories <- split_list(opts$exclude_gene_categories)
if (identical(exclude_categories, "none")) {
  exclude_categories <- character(0)
}

# Pipeline ---------------------------------------------------------------------

# warnings are recorded with the step that raised them; warnings outside a
# timed step come from option checks
warnings_seen <- data.frame(source = character(0), message = character(0))
current_step <- "options"
runtimes <- list()

#' Run one pipeline step, recording its runtime in seconds
timed <- function(step, expr) {
  current_step <<- step
  on.exit(current_step <<- "options")
  start <- Sys.time()
  result <- expr
  runtimes[[step]] <<- round(
    as.numeric(difftime(Sys.time(), start, units = "secs")),
    2
  )
  result
}

run_pipeline <- function() {
  set.seed(opts$seed)
  if (!is.null(opts$malignant_cell_type) && opts$method != "instaprism") {
    warning(
      "--malignant-cell-type only applies to InstaPrism; ignored for ",
      opts$method,
      call. = FALSE
    )
  }
  update_mode <- opts$update_mode %||% "all"
  if (
    update_mode != "pooled" &&
      !is.null(opts$update_mode) &&
      is.null(opts$malignant_cell_type)
  ) {
    warning(
      "--update-mode ",
      update_mode,
      " needs --malignant-cell-type; no reference update runs ",
      "(use --update-mode pooled to update without one)",
      call. = FALSE
    )
  }
  if (update_mode == "pooled" && !is.null(opts$malignant_cell_type)) {
    warning(
      "--update-mode pooled does not use --malignant-cell-type: every ",
      "profile is updated the same way for all bulk samples",
      call. = FALSE
    )
  }

  # reference
  reference <- timed("load_reference", load_reference(opts$reference))
  reference_name <- S4Vectors::metadata(reference)$reference_name %||%
    sub("\\.rds$", "", basename(opts$reference), ignore.case = TRUE)

  # bulk
  bulk <- timed("read_bulk", read_bulk_matrix(opts$bulk))
  message(
    "Read ",
    ncol(bulk$counts),
    " bulk samples; gene IDs detected as ",
    bulk$gene_id_type
  )

  # reference profiles; cell-type-only methods always collapse states
  reference_has_states <- has_cell_states(reference)
  method_states <- method_uses_cell_states(opts$method)
  if (isTRUE(opts$use_cell_states) && !method_states) {
    warning(
      "--use-cell-states ignored: ",
      opts$method,
      " fits cell types only",
      call. = FALSE
    )
  }
  use_cell_states <- method_states &&
    reference_has_states &&
    !isFALSE(opts$use_cell_states)
  profiles <- timed(
    "prepare_profiles",
    prepare_profiles(reference, use_cell_states, opts$min_cells)
  )

  # genes
  harmonised <- timed(
    "harmonise_genes",
    harmonise_genes(profiles, bulk, opts$min_gene_overlap)
  )
  categories <- read_gene_categories(opts$gene_categories)
  harmonised <- timed(
    "exclude_gene_categories",
    exclude_gene_categories(harmonised, categories, exclude_categories)
  )

  method_reference <- list(
    counts = harmonised$reference,
    labels = profiles$labels
  )
  if (method$needs_markers) {
    method_reference$markers <- timed(
      "select_markers",
      select_markers(harmonised$reference, opts$n_markers)
    )
    message("Selected ", length(method_reference$markers), " marker genes")
  }

  # method
  message("Running ", opts$method)
  method_opts <- method_options(opts$method, opts)
  result <- timed(
    "method",
    run_method(opts$method, method_reference, harmonised$bulk, method_opts)
  )

  # outputs
  proportions <- result$proportions
  proportions$reference <- reference_name
  proportions <- proportions[, c(
    "method",
    "reference",
    "level",
    "label",
    "sample_id",
    "quantity",
    "proportion"
  )]

  prefix <- opts$output_prefix %||%
    paste(
      sub("\\.(tsv|csv|rds)$", "", basename(opts$bulk), ignore.case = TRUE),
      reference_name,
      opts$method,
      sep = "_"
    )
  dir.create(opts$output_dir, showWarnings = FALSE, recursive = TRUE)
  out <- function(suffix) file.path(opts$output_dir, paste0(prefix, suffix))

  data.table::fwrite(proportions, out("_proportions.tsv"), sep = "\t")
  if (opts$save_method_objects) {
    saveRDS(result$object, out(paste0("_", opts$method, ".rds")))
  }

  funnel <- c(
    list(bulk_rows_read = bulk$n_genes_read),
    harmonised$funnel
  )
  if (!is.null(result$info$n_genes_used)) {
    funnel$used_by_method <- result$info$n_genes_used
  }

  run_info <- list(
    tool = "deconvolve.R",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    parameters = run_parameters(method_opts),
    method = list(
      name = opts$method,
      options = method_opts,
      quantity = method$quantity,
      levels = unique(proportions$level),
      info = result$info[setdiff(names(result$info), "n_genes_used")]
    ),
    reference = list(
      path = normalizePath(opts$reference),
      name = reference_name,
      format_version = S4Vectors::metadata(
        reference
      )$deconvolution_reference_version,
      build = S4Vectors::metadata(reference)$build,
      has_cell_states = reference_has_states,
      cell_states_used = use_cell_states,
      labels = profiles$labels,
      dropped_labels = profiles$dropped
    ),
    bulk = list(
      path = normalizePath(opts$bulk),
      gene_id_type = bulk$gene_id_type,
      n_samples = ncol(harmonised$bulk),
      samples = I(colnames(harmonised$bulk))
    ),
    genes = list(
      funnel = funnel,
      excluded_categories = I(exclude_categories),
      excluded_genes = lapply(
        split(
          harmonised$excluded_genes$gene,
          harmonised$excluded_genes$category
        ),
        I
      ),
      n_markers = if (method$needs_markers) {
        length(method_reference$markers)
      }
    ),
    batch_dependent_settings = batch_dependent_settings(method_opts),
    warnings = warnings_seen,
    runtimes_seconds = runtimes,
    session_info = session_info_for_json(),
    outputs = list(
      proportions = out("_proportions.tsv"),
      run_info = out("_run_info.json"),
      report = if (!opts$no_report) out("_report.html"),
      method_object = if (opts$save_method_objects) {
        out(paste0("_", opts$method, ".rds"))
      },
      figures = if (!opts$no_report && opts$save_figures) out("_figures")
    )
  )
  jsonlite::write_json(
    hide_paths_in(run_info),
    out("_run_info.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null",
    digits = NA
  )
  message("Wrote ", out("_proportions.tsv"), " and ", out("_run_info.json"))

  if (!opts$no_report) {
    # the proportions and run info are complete by now, so a failed report
    # is a warning, not a failed run
    tryCatch(
      timed(
        "report",
        render_run_report(
          proportions_path = out("_proportions.tsv"),
          run_info_path = out("_run_info.json"),
          output_path = out("_report.html"),
          template = file.path(
            module_dir,
            "templates",
            "run_report_template.Rmd"
          ),
          figures_dir = if (opts$save_figures) out("_figures"),
          figure_format = opts$figure_format
        )
      ),
      error = function(e) {
        warning(
          "Run report failed: ",
          sub("\\.+$", "", conditionMessage(e)),
          ". The proportions and run info are complete; re-render the ",
          "report from them with render_run_report() in R/run_report.R.",
          call. = FALSE
        )
      }
    )
  }
}

#' Command-line parameters for run info, in option order
#'
#' For InstaPrism, `update_mode` shows the mode actually used, including the
#' default when it was not given.
run_parameters <- function(method_opts) {
  parameters <- opts[setdiff(names(opts), "help")]
  if (opts$method == "instaprism") {
    parameters$update_mode <- method_opts$update_mode
  }
  option_order <- vapply(option_list, function(o) o@dest, character(1))
  parameters[order(match(names(parameters), unique(option_order)))]
}

#' Settings in effect that make results depend on the other bulk samples
batch_dependent_settings <- function(method_opts) {
  settings <- stats::setNames(list(), character(0)) # a JSON object when empty
  if (opts$method != "instaprism") {
    return(settings)
  }
  if (isTRUE(method_opts$instaprism_filter)) {
    settings[["instaprism_filter"]] <- paste(
      "InstaPrism's filter drops genes above 1% of reads in more than 10% of",
      "the bulk samples in this run"
    )
  }
  if (
    !is.null(method_opts$malignant_cell_type) &&
      method_opts$update_mode == "all"
  ) {
    settings[["update_mode_all"]] <- paste(
      "The update re-estimates normal cell-type profiles pooled over all",
      "bulk samples in this run"
    )
  }
  if (method_opts$update_mode == "pooled") {
    settings[["update_mode_pooled"]] <- paste(
      "The update re-estimates every cell-type profile pooled over all bulk",
      "samples in this run"
    )
  }
  settings
}

#' Platform and loaded packages from sessioninfo, in a JSON-friendly form
#'
#' Library paths are shown relative to the module directory, and other paths
#' with the home directory replaced by `~`, so run info doesn't expose them.
session_info_for_json <- function() {
  si <- sessioninfo::session_info(info = c("platform", "packages"))
  packages <- as.data.frame(si$packages)
  list(
    platform = lapply(unclass(si$platform), function(x) {
      hide_paths(as.character(x))
    }),
    packages = data.frame(
      package = packages$package,
      version = packages$loadedversion,
      attached = packages$attached,
      date = packages$date,
      source = packages$source,
      library = hide_paths(package_library(packages$package))
    )
  )
}

#' The library (from .libPaths()) each package is installed in
#'
#' renv links packages from its cache, so sessioninfo's own library column is
#' often empty; look the package up in each library instead.
package_library <- function(packages) {
  libs <- .libPaths()
  vapply(
    packages,
    function(pkg) {
      found <- libs[file.exists(file.path(libs, pkg, "DESCRIPTION"))]
      if (length(found) > 0) found[1] else NA_character_
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' Show paths in text relative to the module, and the home directory as `~`
#'
#' Paths under the module directory become relative to it (e.g.
#' `reference-files/gene_categories.tsv`); any other path under the home
#' directory starts with `~`. Works on paths anywhere in the text, such as
#' inside warning messages.
hide_paths <- function(text) {
  text <- gsub(paste0(module_dir, "/"), "", text, fixed = TRUE)
  hide_home(text)
}

#' Apply hide_paths() to every string in a nested list, keeping attributes
#' (so I() vectors stay JSON arrays)
hide_paths_in <- function(x) {
  if (is.list(x)) {
    x[] <- lapply(x, hide_paths_in)
    x
  } else if (is.character(x)) {
    hidden <- hide_paths(x)
    attributes(hidden) <- attributes(x)
    hidden
  } else {
    x
  }
}

#' Replace the home directory in text with `~`
hide_home <- function(text) {
  home <- normalizePath("~", mustWork = FALSE)
  gsub(home, "~", text, fixed = TRUE)
}

invisible(withCallingHandlers(
  run_pipeline(),
  warning = function(w) {
    warnings_seen <<- rbind(
      warnings_seen,
      data.frame(source = current_step, message = conditionMessage(w))
    )
  }
))
