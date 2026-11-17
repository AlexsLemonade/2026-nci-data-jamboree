# Shared code for the held-out evaluation reports:
#   - 05_evaluate-holdout.Rmd (SingleR, scANVI)
#   - 06_compare-bonemarrowmap-holdout.Rmd (BoneMarrowMap vs. SingleR and scANVI)
#
# Source after loading dplyr and ggplot2:
#   source(fs::path(module_dir, "utils", "eval-helpers.R"))
#
# Each report picks its own methods from `method_colors_all`, e.g.
#   method_colors <- method_colors_all[c("SingleR", "scANVI")]
# so a report never gets factor levels or legend entries for methods it doesn't show.

# Plot style -----------------------------------------------------------------------------------

# Color identifies the method; categorical slots 1-3 of the reference palette
method_colors_all <- c(SingleR = "#2a78d6", scANVI = "#eb6834", BoneMarrowMap = "#1baf7a")

ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
grid_color <- "#e1e0d9"
# Single blue ramp for heatmaps
seq_blue <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")

eval_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = ink_primary),
      axis.text = ggplot2::element_text(color = ink_secondary),
      axis.title = ggplot2::element_text(color = ink_secondary),
      strip.text = ggplot2::element_text(color = ink_primary, face = "bold", hjust = 0),
      panel.grid.major = ggplot2::element_line(color = grid_color, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      plot.title.position = "plot"
    )
}

# Metrics --------------------------------------------------------------------------------------

# Per-label precision, recall and F1 within a group of cells. `df` needs `truth` and `pred`
# columns. Labels are those present in the truth; precision counts every cell predicted as the
# label.
per_label_metrics <- function(df) {
  labels <- unique(df$truth)
  purrr::map(labels, \(lab) {
    tp <- sum(df$truth == lab & df$pred == lab)
    n_true <- sum(df$truth == lab)
    n_pred <- sum(df$pred == lab)
    precision <- if (n_pred > 0) tp / n_pred else NA_real_
    recall <- tp / n_true
    f1 <- if (is.na(precision) || precision + recall == 0) 0 else
      2 * precision * recall / (precision + recall)
    tibble::tibble(label = lab, n_cells = n_true, precision, recall, f1)
  }) |>
    purrr::list_rbind()
}

# Accuracy and macro-F1 (unweighted mean of per-label F1) for each method, assay and level.
# `by_level` needs `method`, `assay`, `level` and `correct` columns; `label_metrics` is the
# output of per_label_metrics() grouped by the same columns.
overall_metrics <- function(by_level, label_metrics) {
  by_level |>
    dplyr::group_by(method, assay, level) |>
    dplyr::summarize(n_cells = dplyr::n(), accuracy = mean(correct), .groups = "drop") |>
    dplyr::left_join(
      label_metrics |>
        dplyr::group_by(method, assay, level) |>
        dplyr::summarize(macro_f1 = mean(f1), .groups = "drop"),
      by = c("method", "assay", "level")
    )
}

# Plots ----------------------------------------------------------------------------------------

# Per-label F1 dot plot for one label level, by assay. Point size is the number of true cells.
# `scale_color` is the report's method color scale.
plot_per_label <- function(label_metrics, lvl, scale_color) {
  df <- label_metrics |> dplyr::filter(level == lvl)
  # order labels by mean F1 across methods and assays
  label_order <- df |>
    dplyr::group_by(label) |>
    dplyr::summarize(mean_f1 = mean(f1)) |>
    dplyr::arrange(mean_f1) |>
    dplyr::pull(label)
  df |>
    dplyr::mutate(label = factor(label, levels = label_order)) |>
    ggplot2::ggplot(ggplot2::aes(x = f1, y = label, color = method, size = n_cells)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~assay) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_size_area(max_size = 4, name = "True cells") +
    scale_color +
    ggplot2::labs(x = "F1", y = NULL, title = paste("Per-label F1:", lvl))
}

# Row-normalized confusion heatmap, faceted by assay (rows) and method (columns). `df` needs
# `method`, `assay`, `truth` and `pred` columns. Only shares >= 0.2 are labeled.
plot_confusion <- function(df, title) {
  label_levels <- sort(union(df$truth, df$pred))
  df |>
    dplyr::count(method, assay, truth, pred) |>
    dplyr::group_by(method, assay, truth) |>
    dplyr::mutate(share = n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      truth = factor(truth, levels = rev(label_levels)),
      pred = factor(pred, levels = label_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = pred, y = truth, fill = share)) +
    ggplot2::geom_tile(color = "#fcfcfb", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(share >= 0.2, sprintf("%.2f", share), ""),
        color = share > 0.55
      ),
      size = 2.3,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colors = seq_blue, limits = c(0, 1), name = "Share of true label"
    ) +
    ggplot2::scale_color_manual(values = c(`FALSE` = ink_primary, `TRUE` = "#ffffff")) +
    ggplot2::facet_grid(assay ~ method) +
    ggplot2::labs(x = "Predicted label", y = "True label", title = title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
