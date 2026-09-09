


apptainer exec \
  -B /home -B /omics \
  /home/m168r/decisions/internal/rheinnec/ZP_runtime_improvement.sif \
  Rscript /home/m168r/projects/ZygosityPredictor/inst/tests/full_pid_test.R



apptainer exec \
  -B /home -B /omics \
  /home/m168r/decisions/shared/container/podest_ZP_TAC.sif \
  Rscript /home/m168r/projects/ZygosityPredictor/inst/tests/full_pid_test.R




