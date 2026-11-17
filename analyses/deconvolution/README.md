# Bulk RNA-seq deconvolution

Estimates **cell type** and **cell state** proportions in bulk RNA-seq samples, using a **reference** built from labelled single-cell data.
The core method is [InstaPrism](https://github.com/humengying0907/InstaPrism), with NNLS as a simple alternative.

There are two main scripts:

1. **`build_reference.R`**: labelled single-cell files (SingleCellExperiment `.rds` or AnnData `.h5ad`) → a reference `.rds`.
2. **`deconvolve.R`**: a reference + a bulk count matrix → a proportions table, a run info file and an HTML run report.

Build a reference once, then use it for any number of bulk runs.
The design decisions behind the tool are in [`SPEC.md`](SPEC.md); terms are defined in the repository's [`GLOSSARY.md`](../../GLOSSARY.md).

## Setup

The module uses R 4.6.2 and Bioconductor 3.23, with packages pinned by renv in `renv.lock`.
From this directory, start R (the `.Rprofile` activates renv) and run:

```r
renv::restore()
```

This installs InstaPrism 0.1.6 from GitHub and Biobase, which InstaPrism needs but does not declare.
The run report also needs [pandoc](https://pandoc.org/); without it, the report is skipped with a warning.

## Quick start

```bash
# 1. Build a reference from one or more labelled single-cell files
Rscript build_reference.R \
  --input lib1_processed.rds,lib2_processed.rds \
  --cell-type-column consensus_celltype_annotation \
  --exclude-labels Unknown \
  --output references/my_reference.rds

# 2. Deconvolve bulk samples against it
Rscript deconvolve.R \
  --reference references/my_reference.rds \
  --bulk bulk_counts.tsv \
  --malignant-cell-type Malignant \
  --output-dir results/
```

## Building a reference

`build_reference.R` reads each input file in turn and sums raw counts per **reference sample** × cell state, so memory use stays at about one file.

- **Inputs:** SingleCellExperiment `.rds` files (raw counts in the `counts` assay) or AnnData `.h5ad` files (raw counts in `raw/X`, or in `X` when there is no `raw`). Counts should be raw integers; non-integer values are rounded with a warning. All files must have the same genes in the same order.
- **Labels:** `--cell-type-column` names the cell-type column. `--cell-state-column` optionally adds cell states within cell types; without it, each cell type has a single state of the same name. Cells without a cell type are dropped (`--unlabelled drop`) or labelled `unknown` (`--unlabelled unknown`). `--exclude-labels` drops cells with the listed cell-type or cell-state labels.
- **Reference samples:** taken from `--sample-column` if given, otherwise from a `sample_id` cell column, the sample ID recorded in the file, or the library ID / file name.
- **Output:** a `SummarizedExperiment` of summed counts, one column per reference sample × cell state, plus build provenance and a small cell subsample (up to `--subsample-cap` cells per group, for methods that need single cells).

Small cell types and states are only noted at build time; they are dropped (with a warning) when deconvolving, by `--min-cells`.

## Deconvolving bulk samples

`deconvolve.R` runs one method on one bulk matrix per call. To compare methods or references, run it again and compare the outputs.

**Bulk input:** a genes × samples count matrix as `.tsv` or `.csv` (first column = gene ID) or `.rds` (matrix or data frame with gene IDs as row names).
Gene IDs may be Ensembl IDs or symbols; the type is detected, Ensembl version suffixes are stripped, duplicate IDs are summed and non-integer values are rounded, each with a warning.
Use counts (e.g. salmon estimated counts), not TPM.

**What happens:**

1. Reference profiles are pooled across reference samples, and labels with fewer than `--min-cells` cells are dropped.
2. Bulk genes are matched to the reference. The run warns if fewer than 80% of bulk genes match, and stops below `--min-gene-overlap` (default 0.5).
3. Gene categories from `reference-files/gene_categories.tsv` are excluded: by default mitochondrial, ribosomal and sex-linked genes (`--exclude-gene-categories`; also available: `hemoglobin`, `non_coding`; `none` turns it off).
4. The method runs, and the outputs and report are written.

### Methods

| `--method`             | Levels reported                                         | Notes                                                                     |
| ---------------------- | ------------------------------------------------------- | ------------------------------------------------------------------------- |
| `instaprism` (default) | cell type, and cell state when the reference has states | Optional reference update, below. `--no-cell-states` fits cell types only |
| `nnls`                 | cell type                                               | Non-negative least squares on marker genes (`--n-markers` per cell type)  |

Both estimate each cell type's share of the bulk sample's mRNA (`quantity = mrna_share`), not its share of cells.

### InstaPrism reference update

The single-cell reference may not match the tumour cells in a bulk sample well.
InstaPrism can update the reference to fit the bulk data, and the update is reported as a second fit alongside the initial `instaprism` fit:

| `--update-mode` | Needs `--malignant-cell-type`? | What is updated                                                                | Reported as                   |
| --------------- | ------------------------------ | ------------------------------------------------------------------------------ | ----------------------------- |
| `all` (default) | yes                            | malignant profile per bulk sample; normal profiles once, pooled across the run | `instaprism_update_all`       |
| `malignant`     | yes                            | malignant profile per bulk sample only                                         | `instaprism_update_malignant` |
| `pooled`        | no                             | every profile once, pooled across the run                                      | `instaprism_update_pooled`    |

Without `--malignant-cell-type`, the update runs only in `pooled` mode.

### Results that depend on the whole run

Some settings make one bulk sample's estimates depend on the other samples in the same run:

- **InstaPrism's outlier-gene filter** (on by default; `--no-instaprism-filter` turns it off) drops genes above 1% of reads in more than 10% of the run's bulk samples. With fewer than about 10 samples, a single sample can decide which genes are dropped, and the run warns.
- **`--update-mode all` and `pooled`** re-estimate profiles pooled across the run.

Keep bulk samples that belong together in one run, and don't split a cohort across runs. The report lists the run-dependent settings in effect.

## Outputs

All files go to `--output-dir` and start with `--output-prefix` (default `<bulk file>_<reference name>_<method>`):

| File                       | Contents                                                                                                                                                                                                                            |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `<prefix>_proportions.tsv` | Long table: `method, reference, level, label, sample_id, quantity, proportion`. Proportions sum to 1 per method × level × bulk sample, except `instaprism_update_malignant` cell-state rows, which cover only the normal cell types |
| `<prefix>_run_info.json`   | Parameters, gene counts at each step, excluded and filtered genes, reference summary, warnings, run-dependent settings, runtimes and session info                                                                                   |
| `<prefix>_report.html`     | Self-contained run report (skip with `--no-report`)                                                                                                                                                                                 |
| `<prefix>_<method>.rds`    | The method's raw result, with `--save-method-objects`                                                                                                                                                                               |
| `<prefix>_figures/`        | Each report figure as a file, with `--save-figures` (`--figure-format png` or `pdf`)                                                                                                                                                |

Paths in the run info and report are shown relative to this module, or with the home directory as `~`.

The report is built from the proportions TSV and run info alone, so it can be re-rendered without re-fitting:

```r
source("R/run_report.R")
render_run_report(
  "results/x_proportions.tsv",
  "results/x_run_info.json",
  "results/x_report.html",
  template = "templates/run_report_template.Rmd"
)
```

## Gene category table

`reference-files/gene_categories.tsv` lists genes by category (`gene_id`, `gene_symbol`, `category`), from the Ensembl 104 (GRCh38) annotation that ScPCA uses.
You can supply your own table, with any categories, through `--gene-categories`.

To regenerate the default table:

```bash
bash scripts/download_ensembl_gtf.sh      # to data/genome-references/ at the repository root (git-ignored)
Rscript scripts/make_gene_categories.R
```

## Files

```
build_reference.R                 # build a reference from single-cell files
deconvolve.R                      # deconvolve bulk samples against a reference
R/
  read_inputs.R                   # reading single-cell and bulk inputs
  reference_profiles.R            # building, loading and pooling references
  gene_selection.R                # gene matching, category exclusion, markers
  deconvolution_methods.R         # InstaPrism, NNLS and the method registry
  run_report.R                    # report rendering and figures
templates/run_report_template.Rmd # the run report template
reference-files/gene_categories.tsv
scripts/                          # regenerate the gene category table
SPEC.md                           # design decisions
```
