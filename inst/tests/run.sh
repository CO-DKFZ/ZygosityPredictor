#!/usr/bin/env bash
set -euo pipefail

NEW_IMAGE="${NEW_IMAGE:-/home/m168r/decisions/internal/rheinnec/ZygosityPredictor-devel-1.13.2.sif}"
BASELINE_IMAGE="${BASELINE_IMAGE:-/home/m168r/decisions/shared/container/podest_ZP_TAC.sif}"
PROJECT_DIR="${PROJECT_DIR:-/home/m168r/projects/ZygosityPredictor}"
ZP_WORKERS="${ZP_WORKERS:-1}"

ZP_RUN_MODE="new_workers${ZP_WORKERS}" ZP_WORKERS="${ZP_WORKERS}" apptainer exec \
  -B /home -B /omics \
  "${NEW_IMAGE}" \
  Rscript "${PROJECT_DIR}/inst/tests/full_pid_test.R"

ZP_RUN_MODE="baseline" ZP_WORKERS=1 apptainer exec \
  -B /home -B /omics \
  "${BASELINE_IMAGE}" \
  Rscript "${PROJECT_DIR}/inst/tests/full_pid_test.R"
