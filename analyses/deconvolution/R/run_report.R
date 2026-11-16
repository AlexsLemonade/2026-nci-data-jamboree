# The HTML run report: rendering, data loading and figure functions.
#
# The report is rendered from the proportions TSV and the run info JSON
# alone, so it can be re-rendered without re-fitting:
#   render_run_report("x_proportions.tsv", "x_run_info.json", "x_report.html",
#     template = "templates/run_report_template.Rmd")
# The figure functions below are used by the report template and, with
# --save-figures, write each figure to its own file.

#' Render the run report to a self-contained HTML file
#'
#' @param proportions_path Proportions TSV written by deconvolve.R.
#' @param run_info_path Run info JSON written by deconvolve.R.
#' @param output_path HTML file to write.
#' @param template The report template (`run_report_template.Rmd`).
#' @param figures_dir If given, also save each figure here.
#' @param figure_format `png` or `pdf`, for `figures_dir`.
#' @return `output_path`, invisibly.
render_run_report <- function(
  proportions_path,
  run_info_path,
  output_path,
  template,
  figures_dir = NULL,
  figure_format = "png"
) {
  if (!rmarkdown::pandoc_available()) {
    warning("pandoc not found; skipping the run report", call. = FALSE)
    return(invisible(NULL))
  }
  # render from a scratch copy so knitting leaves nothing in the module
  work_dir <- tempfile("run_report_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  file.copy(template, work_dir)
  if (!is.null(figures_dir)) {
    dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
    figures_dir <- normalizePath(figures_dir)
  }
  output_dir <- normalizePath(dirname(output_path))
  rmarkdown::render(
    file.path(work_dir, basename(template)),
    params = list(
      proportions = normalizePath(proportions_path),
      run_info = normalizePath(run_info_path),
      figures_dir = figures_dir,
      figure_format = figure_format
    ),
    output_file = basename(output_path),
    output_dir = output_dir,
    intermediates_dir = work_dir,
    knit_root_dir = work_dir,
    # the template sees the functions in this file
    envir = new.env(parent = environment(render_run_report)),
    quiet = TRUE
  )
  message("Wrote ", file.path(output_dir, basename(output_path)))
  invisible(output_path)
}

#' Read the proportions TSV and run info JSON for the report
#'
#' @return list(proportions = data.frame, info = list).
load_run_report_data <- function(proportions_path, run_info_path) {
  proportions <- data.table::fread(
    proportions_path,
    sep = "\t",
    data.table = FALSE
  )
  info <- jsonlite::read_json(run_info_path, simplifyVector = TRUE)
  list(proportions = proportions, info = info)
}

# Report helpers ---------------------------------------------------------------

#' The initial and (if any) updated fit names in a proportions table
report_fits <- function(proportions) {
  fits <- unique(proportions$method)
  list(
    initial = setdiff(fits, grep("_update_", fits, value = TRUE))[1],
    updated = grep("_update_", fits, value = TRUE)[1]
  )
}

#' One fixed colour per cell type, the malignant cell type first in red
cell_type_colours <- function(cell_types, malignant = NULL) {
  cell_types <- sort(unique(cell_types))
  has_malignant <- !is.null(malignant) && malignant %in% cell_types
  if (has_malignant) {
    cell_types <- c(malignant, setdiff(cell_types, malignant))
  }
  palette <- c(
    "#4393c3",
    "#5aae61",
    "#fdb863",
    "#8073ac",
    "#e08214",
    "#92c5de",
    "#d6604d",
    "#1b7837",
    "#c51b7d",
    "#35978f",
    "#bababa",
    "#542788"
  )
  n_other <- length(cell_types) - has_malignant
  if (n_other > length(palette)) {
    palette <- scales::hue_pal()(n_other)
  }
  colours <- c(if (has_malignant) "#b2182b", palette[seq_len(n_other)])
  stats::setNames(colours, cell_types)
}

#' Order bulk samples for per-sample plots
#'
#' By decreasing malignant proportion in the initial fit when a malignant
#' cell type is named; otherwise by decreasing proportion of the cell type
#' with the largest mean proportion.
sample_order <- function(proportions, initial_fit, malignant = NULL) {
  ct <- proportions[
    proportions$method == initial_fit & proportions$level == "cell_type",
  ]
  key_type <- if (!is.null(malignant) && malignant %in% ct$label) {
    malignant
  } else {
    means <- tapply(ct$proportion, ct$label, mean, na.rm = TRUE)
    names(means)[which.max(means)]
  }
  key <- ct[ct$label == key_type, ]
  key$sample_id[order(-key$proportion)]
}

#' Format a number with thousands separators
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

#' Show a figure in the report and, with a figures directory, save it
#'
#' @param name File name for the saved figure, without extension.
show_figure <- function(plot, name, figures_dir, figure_format, width, height) {
  print(plot)
  if (!is.null(figures_dir)) {
    ggplot2::ggsave(
      file.path(figures_dir, paste0(name, ".", figure_format)),
      plot,
      width = width,
      height = height,
      dpi = 150
    )
  }
  invisible(plot)
}

# Figures ----------------------------------------------------------------------

x_axis_samples <- function(show_labels) {
  if (show_labels) {
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      panel.grid.major.x = ggplot2::element_blank()
    )
  } else {
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
  }
}

#' Stacked bars of cell-type proportions per bulk sample, for one fit
plot_stacked_proportions <- function(
  ct,
  colours,
  order,
  title,
  subtitle = NULL,
  show_labels = TRUE
) {
  ct$label <- factor(ct$label, levels = names(colours))
  ct$sample_id <- factor(ct$sample_id, levels = order)
  ggplot2::ggplot(ct, ggplot2::aes(sample_id, proportion, fill = label)) +
    ggplot2::geom_col(width = if (show_labels) 0.9 else 1) +
    ggplot2::scale_fill_manual(values = colours, name = "Cell type") +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = if (show_labels) NULL else "Bulk samples",
      y = "Proportion",
      title = title,
      subtitle = subtitle
    ) +
    x_axis_samples(show_labels)
}

#' Cell-state composition within one cell type, per bulk sample
#'
#' @param cs Cell-state rows of one fit for the states of `cell_type`.
plot_cell_states <- function(cs, cell_type, order, show_labels = TRUE) {
  totals <- tapply(cs$proportion, cs$sample_id, sum)
  cs$within <- ifelse(
    totals[cs$sample_id] > 0,
    cs$proportion / totals[cs$sample_id],
    0
  )
  cs$sample_id <- factor(cs$sample_id, levels = order)
  ggplot2::ggplot(cs, ggplot2::aes(sample_id, within, fill = label)) +
    ggplot2::geom_col(width = if (show_labels) 0.9 else 1) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Cell state") +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(
      x = if (show_labels) NULL else "Bulk samples",
      y = "Share within cell type",
      title = paste0("Cell states within ", cell_type)
    ) +
    x_axis_samples(show_labels)
}

#' Heatmap of cell-type proportions (cell types x bulk samples), per fit
plot_proportion_heatmap <- function(ct, colours, order, show_labels = TRUE) {
  ct$label <- factor(ct$label, levels = rev(names(colours)))
  ct$sample_id <- factor(ct$sample_id, levels = order)
  ggplot2::ggplot(ct, ggplot2::aes(sample_id, label, fill = proportion)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~method, ncol = 1) +
    ggplot2::scale_fill_viridis_c(
      labels = scales::percent,
      name = "Proportion"
    ) +
    ggplot2::labs(x = if (show_labels) NULL else "Bulk samples", y = NULL) +
    x_axis_samples(show_labels)
}

#' Updated against initial cell-type proportions, one point per sample x type
plot_fit_comparison <- function(ct, initial_fit, updated_fit, colours) {
  initial <- ct[ct$method == initial_fit, c("sample_id", "label", "proportion")]
  updated <- ct[ct$method == updated_fit, c("sample_id", "label", "proportion")]
  both <- merge(
    initial,
    updated,
    by = c("sample_id", "label"),
    suffixes = c("_initial", "_updated")
  )
  both$label <- factor(both$label, levels = names(colours))
  ggplot2::ggplot(
    both,
    ggplot2::aes(proportion_initial, proportion_updated, colour = label)
  ) +
    ggplot2::geom_abline(linetype = 2, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = colours, name = "Cell type") +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      x = paste0(initial_fit, " (initial fit)"),
      y = updated_fit
    )
}

#' Distribution of cell-type proportions across bulk samples, by fit
plot_fit_distributions <- function(ct, colours) {
  ct$label <- factor(ct$label, levels = names(colours))
  ggplot2::ggplot(ct, ggplot2::aes(label, proportion, fill = method)) +
    ggplot2::geom_boxplot(outlier.size = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_brewer(palette = "Set1", name = "Fit") +
    ggplot2::labs(x = NULL, y = "Proportion")
}

#' Genes remaining at each step, from bulk genes to genes used by the method
plot_gene_funnel <- function(funnel) {
  steps <- c(
    bulk = "bulk genes",
    shared = "shared with reference",
    after_exclusion = "after gene category exclusion",
    used_by_method = "used by the method"
  )
  steps <- steps[names(steps) %in% names(funnel)]
  d <- data.frame(
    step = factor(steps, levels = rev(steps)),
    n = unlist(funnel[names(steps)])
  )
  ggplot2::ggplot(d, ggplot2::aes(n, step)) +
    ggplot2::geom_col(fill = "#4393c3") +
    ggplot2::geom_text(
      ggplot2::aes(label = fmt_n(n)),
      hjust = -0.1,
      size = 3.5
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15)),
      labels = scales::comma
    ) +
    ggplot2::labs(x = "Genes", y = NULL)
}
