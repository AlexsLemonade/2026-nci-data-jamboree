# AML scAtlas reference

Builds a deconvolution **reference** from the AML scAtlas ([Whittle et al., eLife 2026](https://doi.org/10.7554/eLife.104978)), using `analyses/deconvolution/build_reference.R`.
It is a test of reference building on a large, integrated **reference source**, not a vetted AML reference.

## Source data

The [CELLxGENE dataset](https://cellxgene.cziscience.com/collections/071b706a-7ea7-47a4-bddf-6457725839fc) (version `be1f3ced-6b8f-47b0-a3e2-7edc8d466abd`):

- 748,679 cells from 20 studies: 159 AML patients and 51 healthy donors, diagnosis samples only, all bone marrow.
- 21,588 genes, with Ensembl IDs.
- Raw counts in `raw/X`.

The authors also publish their own AnnData object on [figshare](https://doi.org/10.48420/27269946) (CC BY-NC 4.0).
It has more genes (39,075) but uses gene symbols only, and it has no extra malignancy annotation, so it isn't used here.

## Labels

Neither published copy has the authors' cluster-level malignancy call (`malignant_class`), so the reference uses their per-cell calls:

- **`Malignant`** (the **malignant cell type**): cells from AML donors whose `HSPC Cell Type` is a van Galen "-like" label or `LSC`. Those labels are its **cell states**: HSC-like, Prog-like, GMP-like, ProMono-like, Mono-like, cDC-like and LSC.
- **Every other cell** takes its `Author Cell Type` (15 normal-lineage labels), with no cell states.

**Known limitation:** `HSPC Cell Type` is only set for cells whose `Author Cell Type` is `HSPC`.
Leukaemic cells that look like more mature cells keep their normal label.
For example, blasts in monocytic AML are labelled `CD14+ Mono` or `ProMono`.
So the normal myeloid profiles carry leukaemic signal, and the malignant profile under-represents monocytic and more committed AML.
This is acceptable for testing reference building, but don't trust malignant fractions estimated with this reference.

Cells from Smart-seq2 libraries (Velten, Stetson) are dropped, because their counts are reads, not UMIs.
Each `donor_id` is one **reference sample**.

## Usage

Preparation uses Python, with the environment pinned by [uv](https://docs.astral.sh/uv/) (`pyproject.toml`, `uv.lock`).
Run these from this directory:

```bash
# 1. Download the CELLxGENE h5ad to data/aml-scatlas/ at the repository root (git-ignored)
bash scripts/download_aml_scatlas.sh

# 2. Split it into one h5ad per reference sample, with label columns, in data/aml-scatlas/samples/
uv run python scripts/split_by_sample.py

# 3. Build the reference with analyses/deconvolution/build_reference.R, to data/references/aml-scatlas.rds
bash scripts/build_aml_scatlas_reference.sh
```

Step 2 writes `samples.tsv`, which lists each file with its cell and malignant-cell counts.
The output files hold raw counts in `X` and Ensembl IDs with a `gene_symbol` column in `var`.
They add two label columns, `deconv_cell_type` and `deconv_cell_state`.

Step 3 needs the deconvolution module's R environment (see its README).
Building from all 199 reference samples took about 3 minutes, with a peak of about 10 GB of memory.
The result has 21,588 genes, 16 cell types and 22 cell states, covering 744,176 cells.

The source keeps only 21,588 genes, so a typical bulk matrix of about 60k genes shares only about a third of its genes with this reference.
That's fine: `deconvolve.R` measures overlap against the smaller gene set, here the reference's genes, and 99.7% of them are in Ensembl 104.
