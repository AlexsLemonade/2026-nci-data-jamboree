#!/bin/bash
# Download the reference data used by this module:
#   - AML scAtlas from CZI CELLxGENE (SingleR and scANVI references)
#   - BoneMarrowMap Symphony reference and UMAP model (Zeng et al., Blood Cancer Discovery 2025)
#
# Files that already exist are skipped.
#
# Usage (from anywhere):
#   bash analyses/annotate-aml/00_download-reference-data.sh
#
# Environment variables:
#   FORCE_DOWNLOAD  set to 1 to download files again even if they exist (default: 0)

set -euo pipefail

# Run from the directory that this shell script is saved in
cd "$(dirname "$0")"

FORCE_DOWNLOAD=${FORCE_DOWNLOAD:-0}

data_dir="../../data"
scatlas_dir="${data_dir}/AML-scAtlas"
bonemarrowmap_dir="${data_dir}/BoneMarrowMap"
bonemarrowmap_url="https://bonemarrowmap.s3.us-east-2.amazonaws.com"

# Download a URL to a destination file, skipping it if it already exists. Downloads go to a
# temporary file first, so a failed download doesn't leave a partial file that would be skipped
# on the next run.
download() {
  local url="$1"
  local dest="$2"
  if [[ -f "${dest}" && "${FORCE_DOWNLOAD}" != "1" ]]; then
    echo "${dest} exists; skipping (set FORCE_DOWNLOAD=1 to download again)"
    return
  fi
  mkdir -p "$(dirname "${dest}")"
  echo "Downloading ${url} to ${dest}"
  wget -O "${dest}.tmp" "${url}"
  mv "${dest}.tmp" "${dest}"
}

# AML scAtlas from CZI CELLxGENE (~6.4 GB)
download \
  "https://datasets.cellxgene.cziscience.com/be1f3ced-6b8f-47b0-a3e2-7edc8d466abd.h5ad" \
  "${scatlas_dir}/AML_scAtlas.h5ad"

# BoneMarrowMap Symphony reference for query mapping (~344 MB)
download \
  "${bonemarrowmap_url}/BoneMarrowMap_SymphonyReference.rds" \
  "${bonemarrowmap_dir}/BoneMarrowMap_SymphonyReference.rds"

# BoneMarrowMap uwot model, for projecting query cells onto the reference UMAP (~221 MB)
download \
  "${bonemarrowmap_url}/BoneMarrowMap_uwot_model.uwot" \
  "${bonemarrowmap_dir}/BoneMarrowMap_uwot_model.uwot"
