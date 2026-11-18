# Shared code for the inferCNV re-run on SCPCP000007:
#   - 13_run-scpca-infercnv.R (per-gene profiles over unflagged cells)
#   - 14_call-scpca-cnv.R (arm calls and per-cell calls)
#
# See .workspace/annotate-aml/scpcp000007.md (step 4, "Revisiting 14's rules").

# Quantile of the reference cells' noise score above which a cell is flagged as noisy (rule A)
cnv_noise_quantile <- 0.99

# z-scores per cell and arm, relative to the run's reference cells.
#   arm_scores: matrix or data frame, cells x arms (mean residual per arm; no NA)
#   is_reference: logical, one per cell
# Returns a cells x arms matrix with the same column names.
cnv_arm_z <- function(arm_scores, is_reference) {
  arm_scores <- as.matrix(arm_scores)
  stopifnot(!anyNA(arm_scores), length(is_reference) == nrow(arm_scores), any(is_reference))
  mu <- colMeans(arm_scores[is_reference, , drop = FALSE])
  sd <- matrixStats::colSds(arm_scores[is_reference, , drop = FALSE])
  stopifnot("An arm has zero variance in the reference cells" = all(sd > 0))
  z <- sweep(sweep(arm_scores, 2, mu, "-"), 2, sd, "/")
  colnames(z) <- colnames(arm_scores)
  z
}

# Rule A: a cell is noisy if its median |z| over arms is above the cnv_noise_quantile of the
# reference cells. A clonal change affects a few arms and leaves the median |z| low, while a
# shift across most arms in both directions doesn't. Returns a list with the per-cell score,
# the threshold, and the logical flag.
cnv_noise_flag <- function(z, is_reference, quantile = cnv_noise_quantile) {
  score <- matrixStats::rowMedians(abs(z))
  threshold <- stats::quantile(score[is_reference], quantile, names = FALSE)
  list(score = score, threshold = threshold, flagged = score > threshold)
}
