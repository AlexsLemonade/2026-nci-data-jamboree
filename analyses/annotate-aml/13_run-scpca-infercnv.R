#!/usr/bin/env Rscript
# Run inferCNV on the SCPCP000007 AML libraries with a lineage-matched reference from the two
# healthy controls, and summarize each cell's inferred copy-number signal per chromosome arm.
#
# See .workspace/annotate-aml/scpcp000007.md (step 4). ScPCA's own inferCNV used each sample's
# lymphoid cells as the reference, so lineage expression differences between myeloid blasts and
# lymphocytes were read as CNVs (+12q, +1q, ±6p, ... in every library). Here the reference is
# healthy-control cells in five lineage groups. With several reference groups, inferCNV treats
# expression within the range of any normal lineage as copy-number neutral.
#
# Runs:
#   - one per AML library (26): query = the library's cells, reference = both healthy controls
#   - two control runs: each healthy control as the query against a reference from the other
#     only; these should come out flat
#   The AML libraries' own T/NK cells stay in the query as a second check (expected flat).
#
# Cells: singlets (scDblFinder) with >= 500 genes detected (the low-complexity rule). Reference
# cells must also pass project-wide BoneMarrowMap mapping QC. Each reference lineage group is
# downsampled to at most 500 cells.
#
# inferCNV settings: cutoff = 0.1 (10x), denoise = TRUE, HMM = FALSE, no plots. The output
# used is the final residual expression matrix (centered at 1). Per cell and arm:
# arm_<arm> = mean(residual - 1) over the arm's genes (positive = gain-like, negative =
# loss-like), plus cnv_score = mean((residual - 1)^2) over all genes. Calibration against the
# reference cells is done in the report, not here.
#
# Gene positions come from EnsDb.Hsapiens.v86 (GRCh38; the data use Ensembl 104, so a few newer
# genes have no position and are dropped). Autosomes only. Arms are split at approximate GRCh38
# centromere midpoints (UCSC cytoBand "acen" bands, entered by hand, accurate to ~1 Mb). Genes
# within 3 Mb of a centromere are dropped so that error doesn't matter. Arms with fewer than 20
# genes after inferCNV's filtering are NA.
#
# Input:
#   - data/scpca/SCPCP000007/<sample>/<library>_processed.rds
#   - results/scpca/bonemarrowmap/bonemarrowmap_predictions_SCPCP000007.tsv (09)
#   - results/scpca/features/cell_features_SCPCP000007.tsv (10)
#
# Each cell also gets the noise flag from utils/cnv-helpers.R (rule A in the plan): median |z|
# over arms above the 99th percentile of the run's reference cells.
#
# Segment scores, for telling arm-wide shifts from focal ones (rule E): each scored arm is split
# into consecutive segments of about 50 genes (round(n / 50) segments of near-equal size, at
# least 1), in genomic order. That's half of inferCNV's smoothing window (101 genes), so
# neighboring segments overlap in signal. Per cell, a segment's score is the mean residual over
# its genes, the same statistic as the arm scores at finer resolution. (Per-gene medians across
# cells were tried first and dropped: denoising sets most residuals to a constant, so they hide
# shifts that cells carry in different genes.)
#
# Output (results/scpca/infercnv-rerun/):
#   - runs/<run_id>/: inferCNV's working files
#   - arm_scores/<run_id>_arm_scores.tsv: one row per cell in the run, with noise_score and noisy
#   - segment_scores/<run_id>_segment_scores.tsv: one row per cell, one column per segment
#   - segment_scores/<run_id>_segments.tsv: per segment, its arm, index, gene count and span
#   - run_summary.tsv: cells and genes per run
#
# Usage (run from analyses/annotate-aml so renv is activated):
#   Rscript 13_run-scpca-infercnv.R [run_id ...]
# Run IDs are AML library IDs or control_<healthy library ID>. With no arguments, all 28 runs
# are done. Runs with both output files are skipped unless FORCE=1. Threads: N_THREADS (default
# 4).

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
force <- Sys.getenv("FORCE", "0") == "1"
n_threads <- as.integer(Sys.getenv("N_THREADS", "4"))

min_genes <- 500
mad_threshold <- 2.5
max_ref_per_group <- 500
min_genes_per_arm <- 20
segment_genes <- 50
centromere_buffer <- 3e6
seed <- 2026

# Lineage groups (BoneMarrowMap broad labels). Stromal cells aren't used as reference.
lineage_groups <- list(
  HSPC = c("HSC MPP", "LMPP", "Cycling Progenitor", "Early Lymphoid", "MEP", "Early GMP"),
  myeloid = c("Late GMP", "Pro-Monocyte", "Monocyte", "cDC", "pDC", "EoBasoMast Precursor"),
  erythroid_mk = c("Early Erythroid", "Late Erythroid", "Megakaryocyte Precursor"),
  B = c("Pro-B", "Pre-B", "B", "Plasma Cell"),
  T_NK = c("Naive T", "CD4 Memory T", "CD8 Memory T", "NK")
)
lineage_of <- function(broad) {
  out <- rep("other", length(broad))
  for (g in names(lineage_groups)) out[broad %in% lineage_groups[[g]]] <- g
  out
}

# Approximate GRCh38 centromere midpoints (Mb), from UCSC cytoBand acen bands
centromere_mb <- c(
  `1` = 123.4, `2` = 93.9, `3` = 90.9, `4` = 50.0, `5` = 48.8, `6` = 59.8, `7` = 60.1,
  `8` = 45.2, `9` = 43.0, `10` = 39.8, `11` = 53.4, `12` = 35.5, `13` = 17.7, `14` = 17.2,
  `15` = 19.0, `16` = 36.8, `17` = 25.1, `18` = 18.5, `19` = 26.2, `20` = 28.1, `21` = 12.0,
  `22` = 15.0
)

# Paths ----------------------------------------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- fs::path(root_dir, "analyses", "annotate-aml")
project_dir <- fs::path(root_dir, "data", "scpca", "SCPCP000007")
metadata_file <- fs::path(project_dir, "single_cell_metadata.tsv")
scpca_results <- fs::path(module_dir, "results", "scpca")
bmm_file <- fs::path(scpca_results, "bonemarrowmap", "bonemarrowmap_predictions_SCPCP000007.tsv")
features_file <- fs::path(scpca_results, "features", "cell_features_SCPCP000007.tsv")
output_dir <- fs::path(scpca_results, "infercnv-rerun")
runs_dir <- fs::path(output_dir, "runs")
scores_dir <- fs::path(output_dir, "arm_scores")
segments_dir <- fs::path(output_dir, "segment_scores")
fs::dir_create(c(runs_dir, scores_dir, segments_dir))

# cnv_arm_z() and cnv_noise_flag(), shared with 14
source(fs::path(module_dir, "utils", "cnv-helpers.R"))
run_summary_file <- fs::path(output_dir, "run_summary.tsv")

stopifnot(
  "Run 09_predict-scpca-bonemarrowmap.R first" = fs::file_exists(bmm_file),
  "Run 10_extract-scpca-features.R first" = fs::file_exists(features_file)
)

# Libraries and runs ---------------------------------------------------------------------------
libraries <- data.table::fread(metadata_file, data.table = FALSE) |>
  dplyr::filter(diagnosis %in% c("Acute myeloid leukemia", "Non-cancerous")) |>
  dplyr::transmute(
    sample_id = scpca_sample_id,
    library_id = scpca_library_id,
    role = dplyr::if_else(diagnosis == "Non-cancerous", "healthy control", "AML"),
    file = fs::path(project_dir, scpca_sample_id, paste0(scpca_library_id, "_processed.rds"))
  )
healthy_ids <- libraries$library_id[libraries$role == "healthy control"]
stopifnot(length(healthy_ids) == 2, sum(libraries$role == "AML") == 26)

runs <- rbind(
  data.frame(
    run_id = libraries$library_id[libraries$role == "AML"],
    query_library = libraries$library_id[libraries$role == "AML"],
    reference_libraries = paste(healthy_ids, collapse = ",")
  ),
  data.frame(
    run_id = paste0("control_", healthy_ids),
    query_library = healthy_ids,
    reference_libraries = rev(healthy_ids)
  )
)
runs$output_file <- fs::path(scores_dir, paste0(runs$run_id, "_arm_scores.tsv"))
runs$segment_file <- fs::path(segments_dir, paste0(runs$run_id, "_segment_scores.tsv"))
runs$segment_info_file <- fs::path(segments_dir, paste0(runs$run_id, "_segments.tsv"))

if (length(args) > 0) {
  unknown <- setdiff(args, runs$run_id)
  if (length(unknown) > 0) stop("Unknown run IDs: ", paste(unknown, collapse = ", "))
  runs_to_do <- runs[runs$run_id %in% args, ]
} else {
  runs_to_do <- runs
}
if (!force) {
  done <- fs::file_exists(runs_to_do$output_file) & fs::file_exists(runs_to_do$segment_file) &
    fs::file_exists(runs_to_do$segment_info_file)
  runs_to_do <- runs_to_do[!done, ]
}
message(nrow(runs_to_do), " inferCNV runs to do (", n_threads, " threads)")

# Cell selection -------------------------------------------------------------------------------
bmm <- data.table::fread(bmm_file, data.table = FALSE,
                         select = c("cell_id", "library_id", "barcode", "role", "bmm_label_broad",
                                    "bmm_mapping_error"))
features <- data.table::fread(features_file, data.table = FALSE,
                              select = c("cell_id", "detected", "scDblFinder_class"))
qc_threshold <- stats::median(bmm$bmm_mapping_error) + mad_threshold * stats::mad(bmm$bmm_mapping_error)

cells <- bmm |>
  dplyr::inner_join(features, by = "cell_id") |>
  dplyr::mutate(
    lineage = lineage_of(bmm_label_broad),
    eligible = scDblFinder_class == "singlet" & detected >= min_genes,
    qc_pass = bmm_mapping_error < qc_threshold
  )
stopifnot(nrow(cells) == nrow(bmm))

# Reference cells from the given healthy libraries, downsampled per lineage group
select_reference <- function(reference_libraries) {
  set.seed(seed)
  cells |>
    dplyr::filter(library_id %in% reference_libraries, eligible, qc_pass, lineage != "other") |>
    dplyr::group_by(lineage) |>
    dplyr::slice_sample(n = max_ref_per_group) |>
    dplyr::ungroup()
}

# Gene positions and arms ----------------------------------------------------------------------
message("Reading gene positions from EnsDb.Hsapiens.v86")
genes <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                          columns = c("gene_id", "gene_biotype"))
gene_positions <- data.frame(
  gene_id = genes$gene_id,
  chr = as.character(GenomicRanges::seqnames(genes)),
  start = GenomicRanges::start(genes),
  stop = GenomicRanges::end(genes)
) |>
  dplyr::filter(chr %in% names(centromere_mb)) |>
  dplyr::mutate(
    centromere = centromere_mb[chr] * 1e6,
    midpoint = (start + stop) / 2,
    arm = paste0("chr", chr, dplyr::if_else(midpoint < centromere, "p", "q"))
  ) |>
  dplyr::filter(abs(midpoint - centromere) > centromere_buffer) |>
  dplyr::mutate(chr_num = as.integer(chr)) |>
  dplyr::arrange(chr_num, start)
arm_levels <- unique(gene_positions$arm)
message(nrow(gene_positions), " autosomal genes with positions, ", length(arm_levels), " arms")

# Gene order for inferCNV: rownames = genes, columns chr, start, stop
gene_order <- data.frame(
  chr = paste0("chr", gene_positions$chr),
  start = gene_positions$start,
  stop = gene_positions$stop,
  row.names = gene_positions$gene_id
)

# Counts with Ensembl IDs as row names, for the given cells of one library
read_counts <- function(library_id, barcodes) {
  file <- libraries$file[libraries$library_id == library_id]
  sce <- readRDS(file)
  counts <- as(counts(sce), "CsparseMatrix")
  rownames(counts) <- rowData(sce)$gene_ids
  missing <- setdiff(barcodes, colnames(sce))
  if (length(missing) > 0) stop(library_id, ": ", length(missing), " selected barcodes not found")
  counts <- counts[, barcodes, drop = FALSE]
  colnames(counts) <- paste(library_id, barcodes, sep = "-")
  counts
}

# Run ------------------------------------------------------------------------------------------
summaries <- list()

for (i in seq_len(nrow(runs_to_do))) {
  run <- runs_to_do[i, ]
  message("\n[", i, "/", nrow(runs_to_do), "] ", run$run_id)
  reference_libraries <- strsplit(run$reference_libraries, ",")[[1]]

  query_cells <- cells |> dplyr::filter(library_id == run$query_library, eligible)
  reference_cells <- select_reference(reference_libraries)
  # Same order as the count columns below: by reference library, then as selected
  reference_cells <- reference_cells[order(match(reference_cells$library_id, reference_libraries)), ]
  message("  query: ", nrow(query_cells), " cells; reference: ", nrow(reference_cells), " cells (",
          paste(names(table(reference_cells$lineage)), table(reference_cells$lineage), sep = " = ",
                collapse = ", "), ")")

  # Counts: query library plus each reference library, restricted to genes with positions
  counts <- read_counts(run$query_library, query_cells$barcode)
  for (lib in reference_libraries) {
    ref_barcodes <- reference_cells$barcode[reference_cells$library_id == lib]
    counts <- cbind(counts, read_counts(lib, ref_barcodes))
  }
  shared_genes <- intersect(rownames(counts), rownames(gene_order))
  message("  ", length(shared_genes), " of ", nrow(counts), " genes have positions")
  counts <- counts[shared_genes, ]

  annotations <- data.frame(
    group = c(paste0("query_", query_cells$lineage), paste0("reference_", reference_cells$lineage)),
    row.names = c(query_cells$cell_id, reference_cells$cell_id)
  )
  stopifnot(identical(rownames(annotations), colnames(counts)))
  reference_groups <- unique(paste0("reference_", reference_cells$lineage))

  run_dir <- fs::path(runs_dir, run$run_id)
  fs::dir_create(run_dir)

  # Annotation and gene order files in inferCNV's formats (tab-separated, no header)
  annotations_file <- fs::path(run_dir, "cell_annotations.tsv")
  gene_order_file <- fs::path(run_dir, "gene_order.tsv")
  utils::write.table(data.frame(cell = rownames(annotations), group = annotations$group),
                     annotations_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(data.frame(gene = shared_genes, gene_order[shared_genes, ]),
                     gene_order_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annotations_file,
    gene_order_file = gene_order_file,
    ref_group_names = reference_groups
  )
  rm(counts)
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = run_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    no_plot = TRUE,
    no_prelim_plot = TRUE,
    save_rds = FALSE,
    save_final_rds = FALSE,
    num_threads = n_threads
  )

  # Arm scores from the final residual expression
  residual <- as.matrix(infercnv_obj@expr.data) - 1
  gene_arm <- gene_positions$arm[match(rownames(residual), gene_positions$gene_id)]
  stopifnot(!anyNA(gene_arm))
  genes_per_arm <- table(factor(gene_arm, levels = arm_levels))

  arm_scores <- sapply(arm_levels, \(a) {
    rows <- which(gene_arm == a)
    if (length(rows) < min_genes_per_arm) return(rep(NA_real_, ncol(residual)))
    Matrix::colMeans(residual[rows, , drop = FALSE])
  })
  colnames(arm_scores) <- paste0("arm_", arm_levels)

  cell_groups <- annotations[colnames(residual), "group"]
  is_reference <- startsWith(cell_groups, "reference_")

  # Noise flag (rule A) from the arms with scores
  scored_arms <- colnames(arm_scores)[colSums(is.na(arm_scores)) == 0]
  noise <- cnv_noise_flag(cnv_arm_z(arm_scores[, scored_arms, drop = FALSE], is_reference),
                          is_reference)

  scores <- data.frame(
    cell_id = colnames(residual),
    run_id = run$run_id,
    group = cell_groups,
    is_reference = is_reference,
    cnv_score = Matrix::colMeans(residual^2),
    noise_score = noise$score,
    noisy = noise$flagged,
    arm_scores,
    check.names = FALSE,
    row.names = NULL
  )

  # Segment scores (rule E): ~50-gene segments within each scored arm, in genomic order
  gene_start <- gene_positions$start[match(rownames(residual), gene_positions$gene_id)]
  segment_info <- list()
  segment_scores <- list()
  for (a in sub("^arm_", "", scored_arms)) {
    rows <- which(gene_arm == a)
    rows <- rows[order(gene_start[rows])]
    n_segments <- max(1, round(length(rows) / segment_genes))
    segment_of <- ceiling(seq_along(rows) * n_segments / length(rows))
    for (s in seq_len(n_segments)) {
      seg_rows <- rows[segment_of == s]
      name <- sprintf("%s_s%02d", a, s)
      segment_scores[[name]] <- Matrix::colMeans(residual[seg_rows, , drop = FALSE])
      segment_info[[name]] <- data.frame(
        segment = name, arm = a, index = s, n_segments = n_segments, n_genes = length(seg_rows),
        start = min(gene_start[seg_rows]), end = max(gene_start[seg_rows])
      )
    }
  }
  segment_scores <- data.frame(cell_id = colnames(residual), segment_scores, check.names = FALSE,
                               row.names = NULL)
  segment_info <- dplyr::bind_rows(segment_info)

  tmp_file <- paste0(run$output_file, ".tmp")
  data.table::fwrite(scores, tmp_file, sep = "\t", na = "NA")
  fs::file_move(tmp_file, run$output_file)
  tmp_file <- paste0(run$segment_file, ".tmp")
  data.table::fwrite(segment_scores, tmp_file, sep = "\t", na = "NA")
  fs::file_move(tmp_file, run$segment_file)
  data.table::fwrite(segment_info, run$segment_info_file, sep = "\t", na = "NA")
  message("  ", nrow(segment_info), " segments (", min(segment_info$n_genes), "-",
          max(segment_info$n_genes), " genes each)")
  message("  noisy cells: ", sum(noise$flagged), " (",
          paste(names(table(cell_groups[noise$flagged])), table(cell_groups[noise$flagged]),
                sep = " = ", collapse = ", "), ")")
  message("  wrote ", run$output_file, " (", nrow(scores), " cells, ", nrow(residual),
          " genes after filtering; arms with < ", min_genes_per_arm, " genes: ",
          paste(names(genes_per_arm)[genes_per_arm < min_genes_per_arm], collapse = ", "), ")")

  summaries[[run$run_id]] <- data.frame(
    run_id = run$run_id, query_library = run$query_library,
    reference_libraries = run$reference_libraries,
    n_query = nrow(query_cells), n_reference = nrow(reference_cells),
    n_genes_used = nrow(residual)
  )
  rm(infercnv_obj, residual, arm_scores, scores, segment_scores)
  gc(verbose = FALSE)
}

# Summary of the runs done this time, appended to any earlier ones
if (length(summaries) > 0) {
  new_summary <- dplyr::bind_rows(summaries)
  if (fs::file_exists(run_summary_file)) {
    old <- data.table::fread(run_summary_file, data.table = FALSE)
    new_summary <- dplyr::bind_rows(old[!old$run_id %in% new_summary$run_id, ], new_summary)
  }
  data.table::fwrite(new_summary, run_summary_file, sep = "\t")
  message("\nWrote ", run_summary_file)
}

done <- sum(fs::file_exists(runs$output_file) & fs::file_exists(runs$segment_file) &
            fs::file_exists(runs$segment_info_file))
message(done, " of ", nrow(runs), " runs have output")
