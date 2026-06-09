#!/usr/bin/env bash

cd /disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM || exit 1

source /home/lzhenn/.bashrc_liquor_i22wrf415 >/dev/null 2>&1 || true
set -u

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export KMP_INIT_AT_FORK=FALSE

for case_name in india_wb23_pre india_wb23_core india_wb23_post; do
  echo "===== PHASE ${case_name} START $(date)"
  if [[ ! -f "cases/${case_name}/input/base_pres.nc" ||
        ! -f "cases/${case_name}/input/base_surf.nc" ||
        ! -f "cases/${case_name}/input/perturbed_pres.nc" ||
        ! -f "cases/${case_name}/input/perturbed_surf.nc" ]]; then
    python3 -u run_case.py "${case_name}" --step build || exit 10
  else
    echo "Inputs already exist; skipping build for ${case_name}"
  fi
  python3 -u run_case.py "${case_name}" --step run --nproc 160 || exit 20
  python3 -u scripts/validate_additivity.py --case "${case_name}" --plot \
    > "cases/${case_name}/output/validate_additivity.txt" 2>&1 || true
  python3 -u scripts/plot_fig3_self.py --case "${case_name}" || true
  echo "===== PHASE ${case_name} DONE $(date)"
done

python3 -u scripts/plot_fig3_self.py \
  --case india_wb23_pre india_wb23_core india_wb23_post || true
