#!/usr/bin/env Rscript
# Build the default gene category table used by deconvolve.R to exclude gene
# categories before fitting.
#
# Input: an Ensembl GRCh38 GTF (download with download_ensembl_gtf.sh).
# Output: a TSV with one row per gene x category: gene_id, gene_symbol, category.
#
# Categories:
#   mitochondrial  genes on the mitochondrial chromosome (MT-)
#   ribosomal      protein-coding RPS/RPL genes (incl. RPLP0/1/2, RPSA) and
#                  mito-ribosomal MRPS/MRPL
#   sex            all chrY genes, plus XIST and TSIX
#   hemoglobin     HBA1/HBA2/HBB/HBD/HBG1/HBG2/HBE1/HBZ/HBM/HBQ1 (not HBEGF etc.)
#   non_coding     Ensembl gene biotype other than protein_coding
#
# Usage: Rscript make_gene_categories.R [--gtf <path>] [--output <path>]

suppressPackageStartupMessages({
  library(optparse)
})

module_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
})
repo_root <- system2(
  "git",
  c("-C", module_dir, "rev-parse", "--show-toplevel"),
  stdout = TRUE
)

option_list <- list(
  make_option(
    "--gtf",
    default = file.path(
      repo_root,
      "data",
      "references",
      "Homo_sapiens.GRCh38.104.gtf.gz"
    ),
    help = "Ensembl GTF (gzipped) [default: %default]"
  ),
  make_option(
    "--output",
    default = file.path(module_dir, "reference-files", "gene_categories.tsv"),
    help = "Output TSV [default: %default]"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opts$gtf)) {
  stop(
    "GTF not found: ",
    opts$gtf,
    "\nRun scripts/download_ensembl_gtf.sh first."
  )
}

# ---- read gene records from the GTF ----
# Only 'gene' records and the needed attributes are read.
gtf_genes <- Bioc.gff::readGFF(
  opts$gtf,
  columns = "seqid",
  tags = c("gene_id", "gene_name", "gene_biotype"),
  filter = list(type = "gene")
)
genes <- data.frame(
  gene_id = gtf_genes$gene_id,
  gene_symbol = gtf_genes$gene_name,
  biotype = gtf_genes$gene_biotype,
  chrom = as.character(gtf_genes$seqid),
  stringsAsFactors = FALSE
)
stopifnot(!anyDuplicated(genes$gene_id))
message("Read ", nrow(genes), " genes from ", basename(opts$gtf))

# ---- assign categories ----
symbol <- ifelse(is.na(genes$gene_symbol), "", genes$gene_symbol)
category_rules <- list(
  mitochondrial = genes$chrom == "MT",
  # protein-coding ribosomal protein genes only. Names are anchored so that the
  # S6 kinases (RPS6K*), binding proteins (RPS19BP1), readthroughs
  # (RPS10-NUDT3), pseudogenes and antisense transcripts are not included.
  ribosomal = genes$biotype == "protein_coding" &
    grepl(
      "^(RP[SL][0-9]+[A-Z]?[0-9]*|RPL36AL|RPLP[0-9]|RPSA|MRP[SL][0-9]+[A-Z]?)$",
      symbol
    ),
  sex = genes$chrom == "Y" | symbol %in% c("XIST", "TSIX"),
  hemoglobin = grepl("^HB(A1|A2|B|D|G1|G2|E1|Z|M|Q1)$", symbol),
  non_coding = genes$biotype != "protein_coding"
)

categories <- do.call(
  rbind,
  lapply(names(category_rules), function(category) {
    in_category <- category_rules[[category]]
    data.frame(
      gene_id = genes$gene_id[in_category],
      gene_symbol = genes$gene_symbol[in_category],
      category = rep(category, sum(in_category)),
      stringsAsFactors = FALSE
    )
  })
)
categories <- categories[order(categories$category, categories$gene_id), ]

dir.create(dirname(opts$output), showWarnings = FALSE, recursive = TRUE)
write.table(
  categories,
  opts$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

message("Wrote ", nrow(categories), " rows to ", opts$output)
print(table(categories$category))
