#!/bin/bash
#SBATCH -p general
#SBATCH --account=rc_cai_pi
#SBATCH --job-name=scen3_raw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.out
#SBATCH --error=/work/users/a/q/aquijano/mixed_sims/hchs_mixedmodels_simulations/slurm_logs/%x_%A_%a.err

# Usage: sbatch --array=1-50 scripts/slurm/scen3_raw.slurm <population>
# Each task runs BLOCK replicates, cpus-per-task of them at a time (one Stata process per CPU).
set -u
BLOCK=20
N_REPS=1000
POP="$1"

OUT_DIR="${SLURM_SUBMIT_DIR}/data/estimates/${POP}/scen3_raw"
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK + 1 ))
LAST=$(( SLURM_ARRAY_TASK_ID * BLOCK )); [ "$LAST" -gt "$N_REPS" ] && LAST="$N_REPS"

module purge; module add stata/19
export STATATMP="${SLURM_SUBMIT_DIR}/.statatmp"
WORK="${STATATMP}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORK" "$OUT_DIR"
export OUT_DIR WORK POP

run_rep() {
    SIM="$1"
    OUT="${OUT_DIR}/scen3_raw_coef_${SIM}.csv"
    [ -s "$OUT" ] && { echo "rep $SIM skipped"; return; }
    mkdir -p "${WORK}/${SIM}"; cd "${WORK}/${SIM}"
    START=$(date +%s)
    timeout 15m stata-se -b do "${SLURM_SUBMIT_DIR}/scripts/stata/scen3_raw.do" "$SIM" "$POP"
    CODE=$?
    grep -h "^FITCHECK" *.log
    [ -s "$OUT" ] && echo "rep $SIM done in $(( $(date +%s) - START )) s" \
                  || echo "rep $SIM FAILED (exit $CODE) after $(( $(date +%s) - START )) s"
}
export -f run_rep

seq "$FIRST" "$LAST" | xargs -P "$SLURM_CPUS_PER_TASK" -I{} bash -c 'run_rep {}' > "${WORK}/reps.out"
cat "${WORK}/reps.out"
echo "summary: done=$(grep -c ' done in ' "${WORK}/reps.out") skipped=$(grep -c ' skipped' "${WORK}/reps.out") failed=$(grep -c ' FAILED' "${WORK}/reps.out")"

cd "$SLURM_SUBMIT_DIR"; rm -rf "$WORK"
