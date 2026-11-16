# Deconvolution tool: specification

*Status: ready for implementation (first build). Decided 2026-11-16. Terms in **bold** are defined in the repository's [`GLOSSARY.md`](../../GLOSSARY.md).*

## 1. Purpose and scope

This tool estimates **cell type** and **cell state** **proportions** in **bulk samples**. It builds a **reference** from labelled single-cell data and uses InstaPrism as the core method. It is general-purpose: any labelled single-cell data and any bulk count matrix will do. Pediatric AML is the primary use case, but nothing in the tool is AML-specific.

It has two entry points that share functions under `R/`:

1. **`build_reference.R`**: labelled single-cell files → a saved reference (`.rds`).
2. **`deconvolve.R`**: reference + bulk counts → proportions table, run info and HTML report.

**First build:** InstaPrism (with the optional reference update) and NNLS.

**Later additions:** the method interface (§7) is designed to take them.
- MuSiC
- DWLS (solver only, on summed profiles)
- cellGeometry
- dtangle

**Not part of this tool (separate efforts):**
- Creating inputs, i.e. preparing the SCE/h5ad and bulk matrices.
- Choosing the AML reference.
- Benchmarking and evaluation, including the final defaults for `--unlabelled` and `--update-mode`.
- Tests: the first build ships without a test suite.
- Advanced gene-set refinement: the paired bulk/single-cell discordant-gene filter, and a batch-independent alternative to InstaPrism's `filter`.
- Two-step cell-state re-deconvolution of updated malignant profiles.
- References combined from several sources with different labelling schemes.
- MuSiC behaviour with few reference samples, which is deferred until MuSiC is added.

**Prior art:** `claude-science/ccdi-exploration/tall/scripts/31_deconvolve.R` and `30_build_inputs.py`; methods notes in `tall/methods_notes/deconvolution.md`; AML pilot in `pilot/pilot_report.md`. These are mirrored in spirit, not copied: the tool is general, R-only, and has no T-ALL-specific reference ladder.

## 2. Environment

- **R 4.6.2, Bioconductor 3.23**, managed by renv. The pins are in `renv.lock` and `renv/settings.json`.
- **InstaPrism**:
  - Install from GitHub, `humengying0907/InstaPrism`, pinned at v0.1.6, commit `7d3b57cde9`. It is not on CRAN or Bioconductor.
  - Install **Biobase** from Bioconductor first: InstaPrism's `NAMESPACE` imports it, but its `DESCRIPTION` leaves it out.
  - **Attach it with `library(InstaPrism)`.** Calling `InstaPrism::refPrepare()` without attaching fails on `%>%`.
- **Other packages**:
  - SingleCellExperiment, SummarizedExperiment and Matrix;
  - **anndataR** (Bioconductor) and **rhdf5**, installed explicitly since they're only Suggests;
  - nnls, optparse, jsonlite, ggplot2, rmarkdown, knitr.
- **System**: pandoc, for the report.
- **Not used**: zellkonverter, Quarto.

## 3. File layout

```
analyses/deconvolution/
├── build_reference.R
├── deconvolve.R
├── R/
│   ├── read_inputs.R      # SCE/h5ad reading (anndataR + rhdf5), bulk reading, gene ID detection
│   ├── reference.R        # per-library summing, cell subsample, state collapse, --min-cells
│   ├── genes.R            # gene intersection, gene category exclusion, select_markers()
│   ├── methods.R          # run_<method>() functions and the method registry
│   └── report.R           # run report rendering and figure functions
├── templates/run_report.Rmd
├── reference-files/gene_categories.tsv
├── scripts/make_gene_categories.R   # one-off generator for the gene category table
└── README.md
```

Both scripts `source()` the files in `R/`. Outputs go only to a user-chosen directory, never into the module.

## 4. `build_reference.R`

### Options (optparse; lists are comma-separated)

| Option | Default | Meaning |
|---|---|---|
| `--input` | required | One or more `.rds` (SingleCellExperiment) or `.h5ad` files |
| `--cell-type-column` | required | Cell metadata column holding cell-type labels |
| `--cell-state-column` | none | Optional cell metadata column holding cell-state labels |
| `--sample-column` | none | Optional cell metadata column holding reference-sample IDs |
| `--exclude-labels` | none | Cell-type or cell-state labels whose cells are dropped (e.g. `Uncertain`) |
| `--unlabelled` | `drop` | `drop` removes cells with `NA` cell type; `unknown` assigns cell type and cell state `unknown` |
| `--subsample-cap` | `50` | Max cells per reference sample × cell state in the cell subsample; `0` disables it |
| `--seed` | `2026` | Seed for the subsample |
| `--name` | output base name | Reference name, stored in the object and used in outputs |
| `--output` | required | Path of the reference `.rds` |

### Reading inputs

- **SCE `.rds`**:
  - Raw counts are in `counts(sce)`. Ignore the `spliced` and `logcounts` assays.
  - Genes are keyed by Ensembl ID (rownames), with `rowData$gene_ids` and `rowData$gene_symbol`. Many symbols are `NA` or duplicated, so the Ensembl ID is the primary key.
  - In ScPCA files, `sample_id`/`library_id` are **not** in `colData`. They are in `metadata(sce)` (and `metadata(sce)$sample_metadata`), and the reader copies them onto the cells.
- **`.h5ad`**:
  - Read cell metadata (`obs`) and gene metadata with **anndataR**.
  - Read raw counts from `raw/X` with **rhdf5**. In ScPCA processed files, `X` holds logcounts and raw counts exist only in `raw/X` (sparse CSC, cells × genes). Gene metadata is in `raw/var`: index = Ensembl ID, `gene_symbol` is categorical, and code −1 means `NA`.
  - anndataR silently drops `.raw` in every released version as of 2026-11 (scverse/anndataR#259). If a later version supports `.raw`, rhdf5 can be retired.
  - **If the file has no `raw/X`**, read everything, including counts from `X`, with anndataR.
- **In every case, check that counts are integers**, and fail with a clear message otherwise. This catches normalised `X`.

### Reference samples

A cell's reference sample comes from the first of these that exists:
1. `--sample-column`;
2. the sample ID in `metadata(sce)` (or `obs$sample_id`);
3. the input file's base name, with a note in the build log.

Libraries of the same reference sample are summed together.

### Labels

- **Dropping cells:** cells matching `--exclude-labels` are dropped. Cells with an `NA` cell type are handled per `--unlabelled`.
- **No state column:** every cell type has a single cell state named after the cell type.
- **Colliding state names:** InstaPrism needs cell-state names that are unique across cell types. When a state name appears under more than one cell type, those states are stored as `cell_type:cell_state`; all other state names are kept as they are.

### Building the reference

- **One library at a time:** read the library, build a sparse cells × group indicator matrix, sum its counts per (reference sample × cell state) group, record cell counts, draw that library's share of the subsample, then discard it. Peak memory is about one library.
- **Genes:** every input file must have identical gene IDs, or the build fails. Every gene is kept; nothing is filtered at build time.
- **Small labels:** the build log flags cell types and cell states below 10 cells, but nothing is dropped. `--min-cells` is applied in `deconvolve.R`.

## 5. Reference object (`.rds`)

A plain **`SummarizedExperiment`** (no wrapper class):

- **Assay `counts`**: a sparse `dgCMatrix` of summed raw counts, genes × groups.
- **Columns**: one per **reference sample** × cell state. `colData` columns: `sample_id`, `cell_type`, `cell_state`, `n_cells`.
- **Rows**: every input gene. `rowData` columns: `gene_id` (Ensembl, no version suffix) and `gene_symbol`.
- **`metadata(se)`**:
  - `schema_version`: an integer. Start at `1`.
  - `reference_name`
  - `build`: a list holding:
    - input files and their MD5s;
    - the label columns used;
    - the unlabelled mode and the excluded labels;
    - the subsample cap and seed;
    - package versions;
    - a timestamp.
  - `cell_subsample`: a `SingleCellExperiment` with sparse raw counts, the same genes, and `colData` columns `sample_id`, `cell_type`, `cell_state`. It's `NULL` when `--subsample-cap 0`.

Pooled profiles per cell type or cell state are not stored. `deconvolve.R` computes them by summing columns.

`deconvolve.R` accepts any `SummarizedExperiment` with a supported `schema_version` and the required `colData` columns, and refuses anything else with a clear message.

## 6. `deconvolve.R`

### Options

| Option | Default | Meaning |
|---|---|---|
| `--reference` | required | Reference `.rds` from `build_reference.R` |
| `--bulk` | required | Genes × bulk samples count matrix: `.tsv`, `.csv` (first column = gene ID) or `.rds` (matrix/data.frame with gene rownames). The delimiter is taken from the extension |
| `--output-dir` | required | Where outputs are written |
| `--output-prefix` | `<bulk base name>_<reference name>` | Prefix for all output files |
| `--methods` | `instaprism` | Comma-separated: `instaprism`, `nnls` (others once added) |
| `--use-cell-states` / `--no-cell-states` | on if the reference has states | Off collapses cell states to cell types before fitting |
| `--min-cells` | `10` | Cell states/types with fewer cells are dropped, with a warning |
| `--malignant-cell-type` | none | Names the **malignant cell type**; turns on the InstaPrism reference update |
| `--update-mode` | `all` | `all` or `malignant` (see §6, step 6) |
| `--instaprism-filter` / `--no-instaprism-filter` | on | InstaPrism's outlier-gene filter |
| `--exclude-gene-categories` | `mitochondrial,ribosomal,sex` | Categories from the gene category table to exclude; `none` disables |
| `--gene-categories` | `reference-files/gene_categories.tsv` | Gene category table to use |
| `--min-gene-overlap` | `0.5` | Fail if fewer than this fraction of bulk genes match the reference |
| `--n-markers` | `60` | Markers per label for methods that need markers |
| `--n-cores` | `1` | Passed to InstaPrism (`n.core`); does not change results |
| `--seed` | `2026` | |
| `--save-method-objects` | off | Also write each method's raw result object |
| `--no-report` | off | Skip the HTML report |
| `--save-figures` | off | Also write each report figure as a separate file |
| `--figure-format` | `png` | `png` or `pdf` |

Options for one method only are explicit options named with that method's prefix (e.g. `--instaprism-filter`, and later `--cellgeometry-adjust-comp`). There is no generic pass-through. Each `run_<method>()` receives only its own options.

### Pipeline

**1. Load and validate the reference.** See §5.

**2. Read the bulk matrix.**
- **Gene ID type is auto-detected:** if most IDs match `^ENS[A-Z]*G\d+`, they're Ensembl; otherwise symbols.
- **Ensembl version suffixes** (`.12`) are stripped.
- **Duplicate IDs** are summed, with a warning.
- **Non-integer values** (e.g. salmon estimated counts) are rounded, with a warning.

**3. Prepare reference profiles.**
- If cell states are off, collapse each cell type's states into one state named after the type.
- Pool columns across reference samples, giving one summed count profile per cell state (and per cell type).
- Apply `--min-cells` to the pooled labels, after any collapsing. Dropped labels produce a warning at either level, never an error.

**4. Harmonise genes.**
- Match bulk IDs against the reference's `gene_id` (Ensembl) or `gene_symbol`, according to the detected bulk ID type. Duplicate reference symbols are summed, with a warning.
- Intersect reference and bulk genes, and report how many genes were lost.
- Fail if the matched fraction of bulk genes is below `--min-gene-overlap`.
- No other filtering happens here.

**5. Exclude gene categories** (§8). This runs after intersection and before InstaPrism's `filter`, and applies to every method, including marker selection.

**6. Run each requested method** (§7) on the same gene set and bulk matrix.

**InstaPrism:**
- **Build the reference** with the exported `refPrepare(sc_Expr = summed counts per cell state, cell.type.labels, cell.state.labels)`. It accepts summed counts, so per-cell data isn't needed. Don't use internal (`:::`) functions.
- **Initial fit:** `InstaPrism(bulk_Expr, refPhi_cs, filter = <--instaprism-filter>, n.core = <--n-cores>)`.
  - Cell-type proportions are in `@Post.ini.ct@theta`; cell-state proportions are in `@Post.ini.cs@theta`.
  - Report it as method `instaprism`, at level `cell_type`, plus level `cell_state` when cell states are used.
- **Reference update** (only with `--malignant-cell-type`): `InstaPrism_update(obj, bulk_Expr, key = <malignant cell type>, cell.types.to.update = "all"` for `--update-mode all`, or `NULL` for `malignant`).
  - **What it does:** the malignant profile is re-estimated **separately for each bulk sample**. In `all` mode, each normal cell type's profile is also re-estimated **once, pooled over every bulk sample in the initial call**, weighted by sequencing depth.
  - **What it reports:** method `instaprism_updated`.
    - `all` mode: cell-type level only, from `merge_updated_theta()`.
    - `malignant` mode: normal cell types keep their cell-state rows, and the malignant cell type collapses to one row.
  - Updated cell types lose their cell-state proportions, so cell-state proportions for updated types come from the initial fit.
  - **Pass the identical bulk matrix**, same columns in the same order, to both calls. The update matches bulk samples **by column position**, so a reordered or subset matrix runs silently with wrong results.
  - **Look up `@theta` rows by name, not position.** The malignant cell type comes first, and levels are mixed.
- **Never split bulk samples across calls.** `--update-mode all` and `--instaprism-filter` both make results depend on which bulk samples share the run.

**NNLS:**
- Fits at **cell-type level only**, on the markers from `select_markers()`.
- Profiles are the CPM of the summed counts per cell type.
- Each bulk sample, as CPM on the same markers, gets `nnls::nnls()`, and the coefficients are renormalised to sum to 1.
- Quantity is `mrna_share`.

**7. Write outputs** (§9) and render the report unless `--no-report` is set.

## 7. Method interface

- **Signature:** `run_<method>(reference, bulk, options)`. It returns a data frame with columns `method, level, label, sample_id, quantity, proportion`; `deconvolve.R` adds `reference`.
- **Registry:** a method registry in `R/methods.R` declares, for each method:
  - its option names and defaults;
  - the `quantity` it estimates (`mrna_share` or `cell_share`);
  - whether it needs markers or the cell subsample;
  - which levels it reports.
- **Methods without states** (NNLS, and later MuSiC, DWLS, dtangle) always fit at `cell_type` level. Cell-state proportions come only from hierarchy-aware methods (InstaPrism).
- **Notes for later methods:**
  - MuSiC and cellGeometry would use `metadata(se)$cell_subsample`.
  - MuSiC estimates `cell_share`.
  - cellGeometry should report cell-type level as the sum of its states, and expose `adjust_comp`, which is batch-dependent by default.

### `select_markers()`

Uses the claude-science rule. For each label at the fitted level, take the top `--n-markers` genes ranked by log2 CPM in that label minus the maximum log2 CPM over the other labels. Genes must have CPM > 20 and a difference > 0.5. Use the union over labels. It runs after gene category exclusion.

## 8. Gene category exclusion

- **The table:** exclusion is driven by an annotated **gene category table** with columns `gene_id`, `gene_symbol` and `category`. A gene can appear in several rows, one per category.
- **Selecting categories:** `--exclude-gene-categories` picks categories from the table.
- **Your own table:** `--gene-categories` replaces the default table with a user's own, with any categories.
- **Matching:** on Ensembl ID first, falling back to symbol when the bulk is keyed by symbol.

**Default table:** `reference-files/gene_categories.tsv`. The committed script `scripts/make_gene_categories.R` generates it once from Ensembl 104 (GRCh38) annotation, the version ScPCA uses. It's needed because ScPCA gene metadata has no biotype or chromosome.

| Category | Contents | Excluded by default |
|---|---|---|
| `mitochondrial` | `MT-` genes | yes |
| `ribosomal` | RPS/RPL (incl. RP[SL]P-named), RPSA, mito-ribosomal MRPS/MRPL | yes |
| `sex` | all chrY genes, plus XIST and TSIX | yes |
| `hemoglobin` | HBA/HBB/HBD/HBG/HBE/HBZ/HBM/HBQ (not HBEGF) | no: erythroid is a real cell type here |
| `non_coding` | Ensembl biotype other than `protein_coding` | no: opt-in, because it removes ~40k genes |

## 9. Outputs

All files are written to `--output-dir` and start with `--output-prefix`:

| File | Contents |
|---|---|
| `<prefix>_proportions.tsv` | Long table: `method, reference, level, label, sample_id, quantity, proportion`. `level` is `cell_type` or `cell_state`; `quantity` is `mrna_share` or `cell_share`. Proportions sum to 1 per method × level × bulk sample. |
| `<prefix>_run_info.json` | Parameters, detected gene ID type, gene funnel counts and genes removed (per category and by InstaPrism's filter), reference summary (cells/reference samples per label, dropped labels), warnings, bulk sample list, batch-dependent settings in effect, per-step runtimes, package versions. The report is rendered from this file plus the TSV, so it can be re-rendered without re-fitting. |
| `<prefix>_report.html` | Run report (§10); skipped with `--no-report`. |
| `<prefix>_<method>.rds` | Raw method result objects; only with `--save-method-objects`. |
| `<prefix>_figures/` | Each report figure as a separate file, named by its report chunk label (e.g. `proportions_instaprism.png`). Only with `--save-figures`; `--figure-format png` or `pdf`. Made by the same plotting functions as the report. |

## 10. Run report

**Rendering:**
- R Markdown (`templates/run_report.Rmd`), rendered by `deconvolve.R` with `rmarkdown::render()` into **one self-contained HTML file**.
- Use **Quarto-style `#|` chunk options and labels** where possible, to ease a later move to Quarto.
- **Tabbed layout:** methods and plot types are tabs within each section.

**Sections, in order:**

1. **Key facts strip:**
   - number of bulk samples;
   - reference name;
   - number of cell states and cell types;
   - genes used, out of the genes in the bulk matrix;
   - methods;
   - number of warnings.
2. **Warnings**, each with its severity and source:
   - labels dropped by `--min-cells`;
   - rounding;
   - duplicate IDs;
   - low gene overlap;
   - **small bulk sets** (fewer than ~10 bulk samples), where InstaPrism's `filter` lets a single sample decide which genes are dropped.
3. **Batch-dependent settings** active in this run (`--update-mode all`, `--instaprism-filter`), with a one-line reason each and a link to the bulk sample list.
4. **Proportions**:
   - stacked bars per bulk sample for each method, samples ordered by InstaPrism malignant proportion (or by the largest cell type when there's no malignant cell type), one fixed colour per cell type across all plots;
   - cell-state composition within each cell type that has more than one state (InstaPrism initial fit);
   - a heatmap.
   - **Above ~60 bulk samples**, drop per-sample labels and use sorted stacks without labels. Keep the heatmap; per-sample values are in the TSV.
5. **Method comparison**, shown **only when two or more methods ran**:
   - each method against the InstaPrism initial fit;
   - distributions by cell type;
   - a table of mean proportions.
6. **Inputs**:
   - a gene funnel: bulk → shared → after gene category exclusion → after InstaPrism filter;
   - the detected ID type;
   - genes removed per category;
   - a collapsible list of genes removed by InstaPrism's filter;
   - a reference table of cells and reference samples per cell type/state, with dropped labels marked.
7. **Run details**:
   - parameters;
   - a runtime graph (**only when the method-comparison section is shown**);
   - a collapsible bulk sample list;
   - versions;
   - output paths.

**For implementation:** give the cell-state plot a legend per cell type, rather than one shared legend, so it stays readable when there are many states.

A prototype report on synthetic data is on local branch `claude/prototype/run-report` (`analyses/deconvolution/prototype-run-report/`).

## 11. Behaviour notes

**Batch dependence.** Two settings make a bulk sample's estimates depend on the other bulk samples in the run. Both are allowed; the report records them.
- **`--update-mode all`:** the normal cell-type profiles are pooled across the run. In claude-science's T-ALL run, the same reference gave pseudobulk RMSE 0.022 as one 225-sample call versus 0.048 in per-patient calls.
- **InstaPrism's `filter`:** it drops genes above 1% of reads in more than 10% of the run's bulk samples. Estimates shifted by up to about 0.017 between batch compositions.

Only `--update-mode malignant` together with `--no-instaprism-filter` gives estimates for each sample that don't depend on the rest of the run.

**Choosing an update mode.** In the AML pilot:
- `malignant` mode fit high-blast samples better: real-bulk diagnosis RMSE 0.028 versus 0.073 for `all`.
- `malignant` mode overestimated the malignant fraction in low-blast (remission) samples: 0.88 versus 0.33 for `all`, against a truth of about 0.

The default is `all`. The final default is left to benchmarking.

## 12. Decision record

These decisions were made in a planning session on 2026-11-16. The supporting research is on these local branches:
- `claude/research/instaprism-api`: InstaPrism API and update semantics, with source-code citations and synthetic checks;
- `claude/research/scpca-reference-reading`: ScPCA file formats, anndataR, and memory strategy;
- `claude/research/fast-deconvolution-methods`: candidate methods and benchmark evidence.
