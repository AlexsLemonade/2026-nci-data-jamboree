#!/bin/bash
# Predict labels for the held-out study with the held-out SingleR and scANVI references and with
# BoneMarrowMap, then render the evaluation reports.
#
# Prerequisite: the held-out references from run-references.sh
#   RUN_FULL=0 bash analyses/annotate-aml/run-references.sh
#
# Usage (from anywhere):
#   bash analyses/annotate-aml/run-holdout-eval.sh
#
# Environment variables:
#   HOLDOUT_STUDY  held-out study (default: mumme_2023); must match the run-references.sh run
#   N_WORKERS      workers for SingleR classification (default: 4)
#   RUN_SINGLER    set to 0 to skip SingleR predictions (default: 1)
#   RUN_SCANVI     set to 0 to skip scANVI query mapping and predictions (default: 1)
#   RUN_BONEMARROWMAP
#                  set to 0 to skip BoneMarrowMap mapping and predictions (default: 1); skipped
#                  with a message if the BoneMarrowMap reference files haven't been downloaded
#   RUN_REPORT     set to 0 to skip rendering the evaluation reports (default: 1)

set -euo pipefail

# Run from the directory that this shell script is saved in, so renv (.Rprofile) and the
# uv project (pyproject.toml) are picked up
cd "$(dirname "$0")"

HOLDOUT_STUDY=${HOLDOUT_STUDY:-mumme_2023}
N_WORKERS=${N_WORKERS:-4}
RUN_SINGLER=${RUN_SINGLER:-1}
RUN_SCANVI=${RUN_SCANVI:-1}
RUN_BONEMARROWMAP=${RUN_BONEMARROWMAP:-1}
RUN_REPORT=${RUN_REPORT:-1}

subset_file="../../data/AML-scAtlas/AML_scAtlas_le29.h5ad"
singler_ref="results/AML_scAtlas_le29_singler_trained_ref_no-${HOLDOUT_STUDY}.rds"
scanvi_model="results/scanvi_ref_no-${HOLDOUT_STUDY}/scanvi_model"
bmm_ref="../../data/BoneMarrowMap/BoneMarrowMap_SymphonyReference.rds"
bmm_uwot="../../data/BoneMarrowMap/BoneMarrowMap_uwot_model.uwot"
holdout_dir="results/holdout"
report_rmd="05_evaluate-holdout.Rmd"
bmm_report_rmd="06_compare-bonemarrowmap-holdout.Rmd"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

require() {
  if [[ ! -e "$1" ]]; then
    echo "Missing $1; $2" >&2
    exit 1
  fi
}

require "${subset_file}" "run 01_subset-atlas.py (or run-references.sh) first."

# --- SingleR -----------------------------------------------------------------------------------
if [[ "${RUN_SINGLER}" == "1" ]]; then
  require "${singler_ref}" "run RUN_FULL=0 run-references.sh first."
  log "SingleR: classifying ${HOLDOUT_STUDY} cells"
  Rscript 04_predict-holdout-singler.R "${HOLDOUT_STUDY}" "${N_WORKERS}"
fi

# --- scANVI ------------------------------------------------------------------------------------
if [[ "${RUN_SCANVI}" == "1" ]]; then
  require "${scanvi_model}" "run RUN_FULL=0 run-references.sh first."
  log "scANVI: mapping ${HOLDOUT_STUDY} cells (scArches) and predicting"
  uv run python 04_predict-holdout-scanvi.py --holdout-study "${HOLDOUT_STUDY}"
fi

# --- BoneMarrowMap -----------------------------------------------------------------------------
if [[ "${RUN_BONEMARROWMAP}" == "1" ]]; then
  if [[ -f "${bmm_ref}" && -f "${bmm_uwot}" ]]; then
    log "BoneMarrowMap: mapping ${HOLDOUT_STUDY} cells and predicting"
    Rscript 04_predict-holdout-bonemarrowmap.R "${HOLDOUT_STUDY}"
  else
    log "BoneMarrowMap reference files not found; skipping (run 00_download-reference-data.sh)"
  fi
fi

# --- Reports -----------------------------------------------------------------------------------
# Render a report with the held-out study as a parameter
render_report() {
  local rmd="$1"
  local output_file="$2"
  log "Rendering ${rmd}"
  Rscript -e "rmarkdown::render(
    '${rmd}',
    params = list(holdout_study = '${HOLDOUT_STUDY}'),
    output_file = '${output_file}'
  )"
}

if [[ "${RUN_REPORT}" == "1" ]]; then
  render_report "${report_rmd}" "05_evaluate-holdout_${HOLDOUT_STUDY}.nb.html"

  # The BoneMarrowMap report also needs the SingleR and scANVI predictions
  bmm_inputs=(
    "${holdout_dir}/singler_predictions_${HOLDOUT_STUDY}.tsv"
    "${holdout_dir}/scanvi_predictions_${HOLDOUT_STUDY}.tsv"
    "${holdout_dir}/bonemarrowmap_predictions_${HOLDOUT_STUDY}.tsv"
  )
  missing_inputs=0
  for input in "${bmm_inputs[@]}"; do
    [[ -f "${input}" ]] || missing_inputs=1
  done
  if [[ "${missing_inputs}" == "0" ]]; then
    render_report "${bmm_report_rmd}" "06_compare-bonemarrowmap-holdout_${HOLDOUT_STUDY}.nb.html"
  else
    log "Predictions missing for ${bmm_report_rmd}; skipping that report"
  fi
fi

log "Done"
