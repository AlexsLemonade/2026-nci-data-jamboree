#!/bin/bash
# Build the SingleR and scANVI references from the AML scAtlas (donors aged 29 and under).
#
# Each reference is built twice:
#   - full: trained on every study; this is the reference used on ScPCA
#   - held-out: trained without HOLDOUT_STUDY, for validation (outputs suffixed _no-<study>)
#
# Prerequisite: data/AML-scAtlas/AML_scAtlas.h5ad (00_download-reference-data.sh).
#
# Usage (from anywhere):
#   bash analyses/annotate-aml/run-references.sh
#
# Environment variables:
#   HOLDOUT_STUDY  study to leave out for validation (default: mumme_2023)
#   FORCE_SUBSET   set to 1 to rerun 02_subset-atlas.py even if its output exists (default: 0)
#   RUN_FULL       set to 0 to skip the full references (default: 1)
#   RUN_HOLDOUT    set to 0 to skip the held-out references (default: 1)

set -euo pipefail

# Run from the directory that this shell script is saved in, so renv (.Rprofile) and the
# uv project (pyproject.toml) are picked up
cd "$(dirname "$0")"

HOLDOUT_STUDY=${HOLDOUT_STUDY:-mumme_2023}
FORCE_SUBSET=${FORCE_SUBSET:-0}
RUN_FULL=${RUN_FULL:-1}
RUN_HOLDOUT=${RUN_HOLDOUT:-1}

data_dir="../../data/AML-scAtlas"
atlas_file="${data_dir}/AML_scAtlas.h5ad"
subset_file="${data_dir}/AML_scAtlas_le29.h5ad"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Render the SingleR notebook; first argument is the held-out study ("" for none)
render_singler() {
  local holdout="$1"
  local output_file="03_build-singler-ref.nb.html"
  if [[ -n "${holdout}" ]]; then
    output_file="03_build-singler-ref_no-${holdout}.nb.html"
  fi
  Rscript -e "rmarkdown::render(
    '03_build-singler-ref.Rmd',
    params = list(holdout_study = '${holdout}'),
    output_file = '${output_file}'
  )"
}

if [[ ! -f "${atlas_file}" ]]; then
  echo "Missing ${atlas_file}; run 00_download-reference-data.sh first." >&2
  exit 1
fi

# --- Step 1: subset the atlas ------------------------------------------------------------------
if [[ "${FORCE_SUBSET}" == "1" || ! -f "${subset_file}" ]]; then
  log "Step 1: subsetting the atlas"
  uv run python 02_subset-atlas.py
else
  log "Step 1: ${subset_file} exists; skipping (set FORCE_SUBSET=1 to rerun)"
fi

# --- Full references ---------------------------------------------------------------------------
if [[ "${RUN_FULL}" == "1" ]]; then
  log "Step 2: SingleR reference (all studies)"
  render_singler ""

  log "Step 3: scANVI reference (all studies)"
  uv run python 04_train-scanvi-ref.py
fi

# --- Held-out references -----------------------------------------------------------------------
if [[ "${RUN_HOLDOUT}" == "1" ]]; then
  log "Step 2: SingleR reference (without ${HOLDOUT_STUDY})"
  render_singler "${HOLDOUT_STUDY}"

  log "Step 3: scANVI reference (without ${HOLDOUT_STUDY})"
  uv run python 04_train-scanvi-ref.py --holdout-study "${HOLDOUT_STUDY}"
fi

log "Done"
