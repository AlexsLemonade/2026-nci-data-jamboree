#!/bin/bash
# Predict labels for the held-out study with the held-out SingleR and scANVI references, then
# render the evaluation report.
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
#   RUN_REPORT     set to 0 to skip rendering the evaluation report (default: 1)

set -euo pipefail

# Run from the directory that this shell script is saved in, so renv (.Rprofile) and the
# uv project (pyproject.toml) are picked up
cd "$(dirname "$0")"

HOLDOUT_STUDY=${HOLDOUT_STUDY:-mumme_2023}
N_WORKERS=${N_WORKERS:-4}
RUN_SINGLER=${RUN_SINGLER:-1}
RUN_SCANVI=${RUN_SCANVI:-1}
RUN_REPORT=${RUN_REPORT:-1}

subset_file="../../data/AML-scAtlas/AML_scAtlas_le29.h5ad"
singler_ref="results/AML_scAtlas_le29_singler_trained_ref_no-${HOLDOUT_STUDY}.rds"
scanvi_model="results/scanvi_ref_no-${HOLDOUT_STUDY}/scanvi_model"
report_rmd="05_evaluate-holdout.Rmd"

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

# --- Report ------------------------------------------------------------------------------------
if [[ "${RUN_REPORT}" == "1" ]]; then
  if [[ -f "${report_rmd}" ]]; then
    log "Rendering ${report_rmd}"
    Rscript -e "rmarkdown::render(
      '${report_rmd}',
      params = list(holdout_study = '${HOLDOUT_STUDY}'),
      output_file = '05_evaluate-holdout_${HOLDOUT_STUDY}.nb.html'
    )"
  else
    log "${report_rmd} not written yet; skipping the report"
  fi
fi

log "Done"
