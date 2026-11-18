#!/usr/bin/env Rscript
# Call copy-number-altered cells in the SCPCP000007 AML libraries from the lineage-matched inferCNV
# re-run (13_run-scpca-infercnv.R).
#
# See .workspace/annotate-aml/scpcp000007.md (step 4, "Revisiting 14's rules"). The thresholds
# come from the known-normal cells in the data rather than being set by hand:
#
# A. Noisy cells: median |z| over arms above the 99th percentile of the run's reference cells
#    (utils/cnv-helpers.R, shared with 13). They're called "noisy" and left out of the group
#    medians and of the reference cells' score distribution.
# B. Altered arms, per AML library, tested in the HSPC and myeloid groups (>= 20 unflagged
#    cells). All three must hold:
#    1. |median z| is above the largest |median z| in the null pool: every query group of the
#       two control runs, plus the T/NK cells of every AML library (groups with >= 20 unflagged
#       cells).
#    2. |median z - the library's own T/NK median z| is above the largest value of the same
#       difference in the control runs (each non-T/NK group minus T/NK). The library's own
#       normal cells are the baseline, so a shift shared by all of its cells doesn't count.
#    3. The robust z across the 26 AML libraries for that arm and group ((median z - median over
#       libraries) / MAD) is beyond 3.5 (the Iglewicz-Hoaglin outlier convention). A shift seen
#       in most libraries, like the myeloid +19p, reflects cell state rather than copy number.
#    If the HSPC and myeloid groups both pass in opposite directions, the arm isn't called.
# C. Libraries with < 20 unflagged T/NK cells can't run test 2: "no within-library control".
#    The control runs make up the null, so they aren't tested: "control".
# D. Per-cell score: the mean signed z over the library's altered arms. "CNV+" if above the 99th
#    percentile of the run's unflagged reference cells, otherwise "not detected" (not evidence of
#    being normal). Libraries without altered arms: "uninformative". Out-of-sample false-positive
#    checks: the CNV+ rate in the library's own T/NK cells, and in each control run's query
#    cells scored on the same arms and threshold (another donor, another inferCNV run).
#
# Cells left out of inferCNV (doublets, < 500 genes) aren't in the output.
#
# Input (results/scpca/infercnv-rerun/): arm_scores/<run_id>_arm_scores.tsv, run_summary.tsv
# Output (same directory):
#   - cnv_calls_SCPCP000007.tsv: one row per cell per run. Reference cells appear once per run
#     they were used in, so filter on is_reference before joining by cell_id.
#   - cnv_library_summary.tsv: one row per run
#   - cnv_arm_tests.tsv: the three tests per AML library, group (HSPC, myeloid) and arm
#   - cnv_arm_group_summary.tsv: median z over unflagged cells per run, cell group and arm
#   - cnv_thresholds.tsv: the thresholds and the number of values each came from
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 14_call-scpca-cnv.R

suppressPackageStartupMessages({
  library(dplyr)
})

min_group_cells <- 20
call_quantile <- 0.99
outlier_z <- 3.5
tested_groups <- c("query_HSPC", "query_myeloid")
normal_group <- "query_T_NK"

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
cnv_dir <- fs::path(module_dir, "results", "scpca", "infercnv-rerun")
scores_dir <- fs::path(cnv_dir, "arm_scores")
run_summary_file <- fs::path(cnv_dir, "run_summary.tsv")
calls_file <- fs::path(cnv_dir, "cnv_calls_SCPCP000007.tsv")
library_summary_file <- fs::path(cnv_dir, "cnv_library_summary.tsv")
arm_tests_file <- fs::path(cnv_dir, "cnv_arm_tests.tsv")
arm_summary_file <- fs::path(cnv_dir, "cnv_arm_group_summary.tsv")
thresholds_file <- fs::path(cnv_dir, "cnv_thresholds.tsv")

# cnv_arm_z() and cnv_noise_flag(), shared with 13
source(fs::path(module_dir, "utils", "cnv-helpers.R"))

stopifnot("Run 13_run-scpca-infercnv.R first" = fs::file_exists(run_summary_file))
run_summary <- data.table::fread(run_summary_file, data.table = FALSE)
score_files <- fs::path(scores_dir, paste0(run_summary$run_id, "_arm_scores.tsv"))
missing <- run_summary$run_id[!fs::file_exists(score_files)]
if (length(missing) > 0) stop("No arm scores for runs: ", paste(missing, collapse = ", "))
stopifnot(nrow(run_summary) == 28, !anyDuplicated(run_summary$run_id))
run_summary$is_control_run <- startsWith(run_summary$run_id, "control_")
control_runs <- run_summary$run_id[run_summary$is_control_run]
stopifnot(length(control_runs) == 2)

# Pass 1: z-scores, noise flags and group medians per run --------------------------------------
read_run <- function(run_id) {
  scores <- data.table::fread(fs::path(scores_dir, paste0(run_id, "_arm_scores.tsv")),
                              data.table = FALSE)
  stopifnot(all(scores$run_id == run_id), !anyDuplicated(scores$cell_id))

  # Arms that are NA in every cell had too few genes (13); partly NA arms would be a bug
  arm_cols <- grep("^arm_", names(scores), value = TRUE)
  n_na <- colSums(is.na(scores[arm_cols]))
  stopifnot("Arm scores are NA for some cells only" = all(n_na %in% c(0, nrow(scores))))
  arm_cols <- arm_cols[n_na == 0]

  z <- cnv_arm_z(scores[arm_cols], scores$is_reference)
  colnames(z) <- sub("^arm_chr", "", arm_cols)
  noise <- cnv_noise_flag(z, scores$is_reference)
  # 13 computes the same flag from unrounded scores (used for its gene profiles). The TSV has
  # ~15 significant digits, so a cell exactly at the threshold could differ; more is a bug.
  if ("noisy" %in% names(scores)) {
    n_diff <- sum(as.logical(scores$noisy) != noise$flagged)
    if (n_diff > 0) message("  ", run_id, ": noise flag differs from 13 for ", n_diff, " cells")
    stopifnot("Noise flag differs from 13's" = n_diff <= 5)
  }

  group_medians <- lapply(sort(unique(scores$group)), \(g) {
    rows <- scores$group == g & !noise$flagged
    data.frame(
      run_id = run_id,
      group = g,
      is_reference = startsWith(g, "reference_"),
      n_cells = sum(rows),
      n_noisy = sum(scores$group == g & noise$flagged),
      arm = colnames(z),
      median_z = if (sum(rows) > 0) matrixStats::colMedians(z[rows, , drop = FALSE]) else NA_real_,
      row.names = NULL
    )
  }) |> dplyr::bind_rows()

  list(scores = scores, z = z, noise = noise, group_medians = group_medians)
}

runs <- lapply(run_summary$run_id, \(id) {
  message("Reading ", id)
  read_run(id)
})
names(runs) <- run_summary$run_id
group_medians <- dplyr::bind_rows(lapply(runs, `[[`, "group_medians")) |>
  dplyr::mutate(is_control_run = run_id %in% control_runs)

# Pass 2: thresholds from the null cases and the arm tests -------------------------------------
usable <- group_medians |> dplyr::filter(!is_reference, n_cells >= min_group_cells)

# Test 1 null: control-run query groups and AML T/NK groups
null_shift <- usable |> dplyr::filter(is_control_run | group == normal_group)
# Test 2 null: in the control runs, each non-T/NK group minus T/NK on the same arm
tnk_medians <- usable |>
  dplyr::filter(group == normal_group) |>
  dplyr::select(run_id, arm, tnk_median_z = median_z)
null_diff <- usable |>
  dplyr::filter(is_control_run, group != normal_group) |>
  dplyr::inner_join(tnk_medians, by = c("run_id", "arm")) |>
  dplyr::mutate(diff = median_z - tnk_median_z)

shift_threshold <- max(abs(null_shift$median_z))
diff_threshold <- max(abs(null_diff$diff))

thresholds <- data.frame(
  name = c("shift", "shift_vs_T_NK", "cross_library_outlier", "noise_quantile", "call_quantile",
           "min_group_cells"),
  value = c(shift_threshold, diff_threshold, outlier_z, cnv_noise_quantile, call_quantile,
            min_group_cells),
  n_values = c(nrow(null_shift), nrow(null_diff), NA, NA, NA, NA),
  source = c(
    "max |median z|: control-run query groups and AML T/NK groups",
    "max |median z - T/NK median z|: control-run non-T/NK groups",
    "Iglewicz-Hoaglin convention for robust z",
    "reference-cell quantile of median |z| over arms (noise flag)",
    "reference-cell quantile of the altered-arm score (CNV+ threshold)",
    "minimum unflagged cells for a group median"
  )
)
message("\nThresholds:")
print(thresholds[, c("name", "value", "n_values")], row.names = FALSE)

# Tests for the AML libraries' HSPC and myeloid groups
aml_tested <- usable |>
  dplyr::filter(!is_control_run, group %in% tested_groups)
cross_library <- aml_tested |>
  dplyr::group_by(group, arm) |>
  dplyr::summarize(n_libraries = dplyr::n(), library_median = stats::median(median_z),
                   library_mad = stats::mad(median_z), .groups = "drop")
stopifnot("Cross-library MAD is zero" = all(cross_library$library_mad > 0))

arm_tests <- aml_tested |>
  dplyr::left_join(tnk_medians, by = c("run_id", "arm")) |>
  dplyr::left_join(cross_library, by = c("group", "arm")) |>
  dplyr::mutate(
    diff_vs_T_NK = median_z - tnk_median_z,
    cross_library_z = (median_z - library_median) / library_mad,
    pass_shift = abs(median_z) > shift_threshold,
    pass_vs_T_NK = !is.na(diff_vs_T_NK) & abs(diff_vs_T_NK) > diff_threshold,
    pass_outlier = abs(cross_library_z) > outlier_z,
    pass_all = pass_shift & pass_vs_T_NK & pass_outlier,
    # Direction must agree with the T/NK-relative and cross-library shifts
    pass_all = dplyr::coalesce(pass_all & sign(diff_vs_T_NK) == sign(median_z) &
                                 sign(cross_library_z) == sign(median_z), FALSE)
  ) |>
  dplyr::select(run_id, group, arm, n_cells, n_noisy, median_z, tnk_median_z, diff_vs_T_NK,
                library_median, library_mad, cross_library_z, pass_shift, pass_vs_T_NK,
                pass_outlier, pass_all)

# Altered arms per library; conflicting directions between groups aren't called
altered_by_run <- arm_tests |>
  dplyr::filter(pass_all) |>
  dplyr::group_by(run_id, arm) |>
  dplyr::summarize(
    direction = if (dplyr::n_distinct(sign(median_z)) == 1) sign(median_z[1]) else 0,
    groups = paste(sub("^query_", "", group), collapse = ","),
    .groups = "drop"
  )
conflicts <- altered_by_run |> dplyr::filter(direction == 0)
if (nrow(conflicts) > 0) {
  message("\nArms with opposite directions in HSPC and myeloid (not called):")
  print(conflicts, row.names = FALSE)
}
altered_by_run <- altered_by_run |> dplyr::filter(direction != 0)

arm_order <- colnames(runs[[1]]$z)
altered_arms_of <- function(run_id) {
  a <- altered_by_run[altered_by_run$run_id == run_id, ]
  a <- a[order(match(a$arm, arm_order)), ]
  stats::setNames(a$direction, a$arm)
}
format_arms <- function(direction) {
  if (length(direction) == 0) return("")
  paste0(ifelse(direction > 0, "+", "-"), names(direction), collapse = ";")
}

# Pass 3: per-cell calls -----------------------------------------------------------------------
arm_score <- function(z, direction) {
  as.vector(z[, names(direction), drop = FALSE] %*% direction) / length(direction)
}

call_run <- function(run) {
  r <- runs[[run$run_id]]
  scores <- r$scores
  ref <- scores$is_reference
  noisy <- r$noise$flagged
  n_unflagged <- function(g) sum(scores$group == g & !noisy)
  n_tnk <- n_unflagged(normal_group)
  direction <- if (run$is_control_run) numeric(0) else altered_arms_of(run$run_id)

  status <- dplyr::case_when(
    run$is_control_run ~ "control",
    n_tnk < min_group_cells ~ "no within-library control",
    length(direction) > 0 ~ "informative",
    TRUE ~ "no altered arms"
  )

  fpr_controls <- stats::setNames(rep(NA_real_, length(control_runs)),
                                  paste0("fpr_", control_runs))
  if (length(direction) > 0) {
    cell_score <- arm_score(r$z, direction)
    threshold <- stats::quantile(cell_score[ref & !noisy], call_quantile, names = FALSE)
    call <- dplyr::case_when(
      ref ~ "reference",
      noisy ~ "noisy",
      cell_score > threshold ~ "CNV+",
      TRUE ~ "not detected"
    )
    # Out-of-sample check on the control runs' unflagged query cells
    for (cr in control_runs) {
      c_run <- runs[[cr]]
      keep <- !c_run$scores$is_reference & !c_run$noise$flagged
      fpr_controls[paste0("fpr_", cr)] <-
        mean(arm_score(c_run$z[keep, , drop = FALSE], direction) > threshold)
    }
  } else {
    cell_score <- NA_real_
    threshold <- NA_real_
    call <- dplyr::case_when(ref ~ "reference", noisy ~ "noisy", TRUE ~ "uninformative")
  }

  calls <- data.frame(
    run_id = run$run_id,
    cell_id = scores$cell_id,
    library_id = sub("-.*$", "", scores$cell_id),
    group = scores$group,
    is_reference = ref,
    cnv_score = scores$cnv_score,
    noise_score = r$noise$score,
    noisy = noisy,
    cnv_arm_score = cell_score,
    cnv_threshold = threshold,
    cnv_call = call
  )

  frac_noisy <- function(g) {
    rows <- scores$group == g
    if (sum(rows) == 0) return(NA_real_)
    mean(noisy[rows])
  }
  frac_cnv <- function(g) {
    rows <- calls$group == g & !noisy
    if (sum(rows) == 0 || length(direction) == 0) return(NA_real_)
    mean(calls$cnv_call[rows] == "CNV+")
  }
  library_summary <- data.frame(
    run_id = run$run_id,
    query_library = run$query_library,
    is_control_run = run$is_control_run,
    n_query = sum(!ref),
    n_reference = sum(ref),
    n_HSPC = n_unflagged("query_HSPC"),
    n_myeloid = n_unflagged("query_myeloid"),
    n_T_NK = n_tnk,
    frac_noisy_query = mean(noisy[!ref]),
    frac_noisy_HSPC = frac_noisy("query_HSPC"),
    frac_noisy_myeloid = frac_noisy("query_myeloid"),
    frac_noisy_T_NK = frac_noisy("query_T_NK"),
    noise_threshold = r$noise$threshold,
    library_status = status,
    altered_arms = format_arms(direction),
    cnv_threshold = threshold,
    frac_cnv_HSPC = frac_cnv("query_HSPC"),
    frac_cnv_myeloid = frac_cnv("query_myeloid"),
    frac_cnv_erythroid_mk = frac_cnv("query_erythroid_mk"),
    frac_cnv_B = frac_cnv("query_B"),
    frac_cnv_T_NK = frac_cnv("query_T_NK"),
    t(fpr_controls),
    check.names = FALSE
  )
  list(calls = calls, library_summary = library_summary)
}

results <- lapply(seq_len(nrow(run_summary)), \(i) call_run(run_summary[i, ]))
calls <- dplyr::bind_rows(lapply(results, `[[`, "calls"))
library_summary <- dplyr::bind_rows(lapply(results, `[[`, "library_summary"))

# Each query cell appears in exactly one run
query_ids <- calls$cell_id[!calls$is_reference]
stopifnot("Query cells appear in more than one run" = !anyDuplicated(query_ids))

arm_summary <- group_medians |>
  dplyr::left_join(altered_by_run |> dplyr::transmute(run_id, arm, altered = TRUE),
                   by = c("run_id", "arm")) |>
  dplyr::mutate(altered = dplyr::coalesce(altered, FALSE))

data.table::fwrite(calls, calls_file, sep = "\t", na = "NA")
data.table::fwrite(library_summary, library_summary_file, sep = "\t", na = "NA")
data.table::fwrite(arm_tests, arm_tests_file, sep = "\t", na = "NA")
data.table::fwrite(arm_summary, arm_summary_file, sep = "\t", na = "NA")
data.table::fwrite(thresholds, thresholds_file, sep = "\t", na = "NA")
message("\nWrote ", calls_file, " (", nrow(calls), " rows, ", length(query_ids), " query cells)")
message("Wrote ", library_summary_file, ", ", arm_tests_file, ", ", arm_summary_file, ", ",
        thresholds_file)

# Summary --------------------------------------------------------------------------------------
message("\nLibrary status:")
print(table(library_summary$library_status))
message("\nTests passed (AML HSPC/myeloid group x arm tests: ", nrow(arm_tests), "):")
print(colSums(arm_tests[c("pass_shift", "pass_vs_T_NK", "pass_outlier", "pass_all")]))
message("\nInformative libraries:")
library_summary |>
  dplyr::filter(library_status == "informative") |>
  dplyr::select(run_id, altered_arms, frac_noisy_HSPC, frac_cnv_HSPC, frac_cnv_myeloid,
                frac_cnv_T_NK, dplyr::starts_with("fpr_")) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with(c("frac", "fpr")), \(x) round(x, 3))) |>
  print(row.names = FALSE)
message("\nQuery cells by call:")
print(table(calls$cnv_call[!calls$is_reference]))
