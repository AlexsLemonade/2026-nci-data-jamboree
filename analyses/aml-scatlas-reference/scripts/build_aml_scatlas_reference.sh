#!/bin/bash
# Build the AML scAtlas deconvolution reference from the per-sample h5ad files
# written by split_by_sample.py, using analyses/deconvolution/build_reference.R.
#
# Writes data/references/aml-scatlas.rds at the repository root (git-ignored).
#
# Usage: bash build_aml_scatlas_reference.sh

set -euo pipefail

repo_root="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
samples_dir="${repo_root}/data/aml-scatlas/samples"
out_dir="${repo_root}/data/references"
mkdir -p "${out_dir}"

# one file per reference sample, as listed in the split manifest
inputs="$(tail -n +2 "${samples_dir}/samples.tsv" | cut -f2 \
  | sed "s|^|${samples_dir}/|" | paste -sd, -)"

# run from the deconvolution module so its renv library is used
cd "${repo_root}/analyses/deconvolution"
Rscript build_reference.R \
  --input "${inputs}" \
  --cell-type-column deconv_cell_type \
  --cell-state-column deconv_cell_state \
  --sample-column donor_id \
  --output "${out_dir}/aml-scatlas.rds"
