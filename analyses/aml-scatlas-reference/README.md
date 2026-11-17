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

## Usage

```bash
bash scripts/download_aml_scatlas.sh   # to data/aml-scatlas/ at the repository root (git-ignored)
```
