#!/bin/bash
# Download the AML scAtlas (Whittle et al., eLife 2026,
# doi:10.7554/eLife.104978) from CELLxGENE: one h5ad with raw counts in raw/X
# and Ensembl gene IDs. The download is pinned to a dataset version and checked
# against its published byte size.
#
# Writes data/aml-scatlas/aml-scatlas_cellxgene.h5ad at the repository root
# (git-ignored). Skipped if the file is already present.
#
# Usage: bash download_aml_scatlas.sh

set -euo pipefail

# CELLxGENE collection 071b706a-7ea7-47a4-bddf-6457725839fc,
# dataset 4568ebb2-468c-4381-9ce7-81cd97abcd22
version="be1f3ced-6b8f-47b0-a3e2-7edc8d466abd"
url="https://datasets.cellxgene.cziscience.com/${version}.h5ad"
expected_size=6398033807
file="aml-scatlas_cellxgene.h5ad"

repo_root="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
out_dir="${repo_root}/data/aml-scatlas"
mkdir -p "${out_dir}"
path="${out_dir}/${file}"

if [[ -f "${path}" ]]; then
  echo "Already present: ${path}"
  exit 0
fi

echo "Downloading ${url}"
curl -fL --retry 3 -o "${path}.tmp" "${url}"
observed="$(wc -c < "${path}.tmp" | tr -d ' ')"
if [[ "${observed}" != "${expected_size}" ]]; then
  echo "Size mismatch for ${file}: expected ${expected_size}, got ${observed}" >&2
  exit 1
fi
mv "${path}.tmp" "${path}"
echo "Size OK (${observed} bytes): ${path}"
