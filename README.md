# hchs_mixedmodels_simulations

Simulation study of survey-weighted three-level linear mixed models (visits in
subjects in households) under the HCHS/SOL sampling design.

## Model

```
y_gfr ~ age_strat_new + x12 + x17 + x18 + x6 + y_bmi + (1 | household) + (1 | subject)
```

Targets: the fixed effects, the household and subject variances (`var_hhid`,
`var_subid`) and the residual variance (`s2`). The truth is the same model fitted,
unweighted, to the full population.

## Methods

All scenarios use robust (sandwich) SEs clustered on households.

**Weight scaling** (Pfeffermann et al., 1998). For a weight `w_ij` of unit *i*
in cluster *j* of the level above:

- `s1`, effective sample size: `w_ij × Σw_ij / Σw_ij²`. The scaled weights sum to the effective cluster size.
- `s2`, cluster sample size: `w_ij × n_j / Σw_ij`. The scaled weights sum to the actual cluster size.

**scen0: one-step pseudo-likelihood, global scaling** (complete data)
The household keeps its unscaled weight `w_bghh`. The subject weight is scaled
globally to mean 1 across the sample, with no scaling within household.

**scen1 / scen2: one-step pseudo-likelihood** (complete data)
Each level except the upper one is scaled with respect to the level above. The
household (upper level) keeps its unscaled weight `w_bghh`. The subject weight is
scaled within household, with `s1` in scen1 and `s2` in scen2.

**scen3: two-step pseudo-likelihood** (complete data)
1. Fit with globally scaled subject weights to estimate the fixed effects and `s2`.
2. Hold those fixed and refit with the household weight and the `s1`-scaled subject
   weight to estimate `var_hhid` and `var_subid`.

Why two steps: in a single fit, global scaling (scen0) leaves the fixed effects
unbiased but biases the variance components. Scaling within household (scen1,
scen2) fixes the variance components but biases the fixed effects. scen3 takes
each part from the scaling that estimates it well.

**scen4 / scen5: multiple imputation + scen3**
Visit 2–3 outcomes are set missing at random and imputed by chained regression
(M = 5), with log design weights as imputation covariates. Each completed dataset
is fitted with scen3, and the results are pooled with Rubin's rules. scen4 uses
`s1`-scaled weights as covariates and scen5 uses `s2`-scaled weights:

| Variant | Imputation weight covariates |
|---|---|
| v1 | `log(w_bg) + log(hh_s)` |
| v2 | `log(w_bg) + log(hh_s) + log(sub_s)` |
| v3 | `log(w_bg) + log(hhsub_s)` |

**Weights by run, step and level.** Visits always have weight 1. `w_bghh` is the raw
household weight. `w_sub_glob` is the subject weight `W_sub` scaled to mean 1 over the
sample. `sub_s1` and `sub_s2` are `W_sub` scaled within household. The fixed effects and
`s2` come from step 1, and `var_hhid` and `var_subid` from step 2. A one-step run takes
every estimate from its single fit.

| Run | Step 1 household | Step 1 subject | Step 2 household | Step 2 subject |
|---|---|---|---|---|
| scen0 | `w_bghh` | `w_sub_glob` | — (one step) | — |
| scen1 | `w_bghh` | `sub_s1` | — (one step) | — |
| scen2 | `w_bghh` | `sub_s2` | — (one step) | — |
| scen3, scen3_strata, scen4/5 (each imputation) | `w_bghh` | `w_sub_glob` | `w_bghh` | `sub_s1` |
| scen3_s2 | `w_bghh` | `w_sub_glob` | `w_bghh` | `sub_s2` |
| scen3_vunw `hh` | `w_bghh` | `w_sub_glob` | `w_bghh` | 1 |
| scen3_vunw `none` | `w_bghh` | `w_sub_glob` | 1 | 1 |
| scen3_raw | `w_bghh` | `w_sub_glob` | `w_bghh` | `w_sub` (raw) |
| scen3_unw | 1 | 1 | — (one step) | — |
| scen3_srs (SRS of households) | 1 | 1 | — (one step) | — |

**Diagnostics for scen3 on oct2026** (complete data)
On oct2026, scen3 overestimates `var_subid` (+3.4%, coverage 0.83), while on dec2024 it is
unbiased. Five variants locate the cause:

| Run | Change from scen3 | `var_subid` bias / coverage |
|---|---|---|
| scen3_strata | Design-stratum indicators in the model and in the truth | +3.4% / 0.83 |
| scen3_srs | Simple random sample of households from the population, unweighted | 0.0% / 0.95 |
| scen3_unw | Same design samples as scen3, unweighted | +0.3% / 0.94 |
| scen3_vunw `none` | Step 1 as scen3; step 2 unweighted at both levels | +0.4% / 0.94 |
| scen3_vunw `hh` | Step 1 as scen3; step 2 with household weight `w_bghh` only, subjects unweighted | +3.1% / 0.85 |

The bias comes from the step 2 weights, not from the sample selection or the population.
The household weight alone reproduces almost all of it (scen3_vunw `hh`), so the subject
weight in step 2 adds little. The fixed effects still need the weights: unweighted,
`age_strat_new` is 17% off (coverage 0.24). With step 1 weighted and step 2 unweighted
(scen3_vunw `none`), the fixed effects stay unbiased (coverage 0.94–0.96).

**Why oct2026 and not dec2024: informative subject weights.** The two populations'
samples have the same design: 1.37 subjects per household, 69.5% single-subject households,
the same within-household variation of `W_sub` and the same `s1`/`s2` scaled weights. They differ
in how informative the subject weight is. Within households, the correlation between log `W_sub` and the
subject's deviation from the household mean outcome is 0.40 on oct2026 against 0.21 on dec2024
(`scripts/evaluation/compare_weights.R`, 50 samples per population,
`data/summaries/diagnostics/weight_comparison_dec2024_oct2026.csv`). This uses the raw outcome,
not residuals, so part of it may be age entering the model only through `age_strat_new`.

**Subject-weight scaling on oct2026.** scen3_s2 is scen3 with the `s2`-scaled subject weight in step 2.

| Run | Subject weight | `var_hhid` bias / coverage | `var_subid` bias / coverage | Fixed effects |
|---|---|---|---|---|
| scen3 | `sub_s1`, step 2 | −7.8% / 0.97 | +3.4% / 0.83 | unbiased |
| scen3_s2 | `sub_s2`, step 2 | +9.9% / 0.88 | +1.1% / 0.93 | unbiased |
| scen1 | `sub_s1`, one step | −7.8% / 0.97 | +3.1% / 0.86 | biased (coverage 0.89–0.91) |
| scen2 | `sub_s2`, one step | +9.9% / 0.88 | +0.8% / 0.93 | biased (coverage 0.88–0.91) |
| scen3_unw | none | −7.1% / 0.96 | +0.3% / 0.94 | `age_strat_new` 17% off |
| scen3_vunw `none` | none in step 2 (household unweighted too) | −6.2% / 0.96 | +0.4% / 0.94 | unbiased |
| scen3_vunw `hh` | none in step 2 (household weighted) | +1.2% / 0.94 | +3.1% / 0.85 | unbiased |
| scen3_raw | raw `w_sub`, step 2 | +624% / 0.00 | −72% / 0.00 | unbiased |

`s1` and `s2` give the same `var_hhid + var_subid` (202.1 against a truth of 197.9, +2.1%).
The weighting inflates the total between-subject variance, and the scaling only decides how the
excess is split between the household and subject levels. Neither scaling removes it. Dropping
the subject weight from step 2 does not remove it either: with the household weight alone the
total is 203.7 (+2.9%), and with no weights in step 2 it is 197.3 (−0.3%).

With the raw subject weight in step 2 (scen3_raw), most of the subject variance moves to the
household level (`var_hhid` 164.6 against 22.7, `var_subid` 49.6 against 175.2), and the total
rises to 214.2 (+8.2%). Scaling the subject weight within household is necessary; the remaining
`var_subid` bias in scen3 comes from the household weight.

## Runs

| Scenario | Populations |
|---|---|
| scen0 | dec2024 |
| scen1, scen2 | dec2024, oct2026 |
| scen3 | dec2024, oct2026 |
| scen3_strata | dec2024, oct2026 |
| scen3_srs, scen3_unw, scen3_s2 | oct2026 |
| scen3_vunw `none`, `hh` | oct2026 |
| scen3_raw | oct2026 |
| scen4 v1–v3, scen5 v1–v3 | oct2026 |

## Layout

```
scripts/
  stata/        scen{0,1,2,3}.do, scen3_{strata,srs,unw,s2,vunw,raw}.do, scen{4,5}_mi.do,
                population_truth_hhsub.do, population_truth_strata.do
  slurm/        population_truth.slurm, scen{0,1,2,3}.slurm, scen3_{strata,srs,unw,s2,vunw,raw}.slurm,
                scen{4,5}_mi.slurm, submit_coverage.slurm
  evaluation/   coverage_assessment.R, compare_weights.R
data/
  population_coefficients/   truth per population and model: _hhsub, _strata (tracked)
  estimates/<pop>/<scen>/    one CSV per replicate (not tracked)
  summaries/<pop>/           coverage tables (tracked)
  summaries/diagnostics/     weight comparison across populations (tracked)
```

## Running

Submit everything from the repo root. Logs go to `slurm_logs/`.

**1. Population truth**, one job per population:

```
sbatch scripts/slurm/population_truth.slurm dec2024
sbatch scripts/slurm/population_truth.slurm oct2026
sbatch scripts/slurm/population_truth.slurm oct2026 strata    # truth for scen3_strata
```

**2. Time one replicate.** Submit only the first array task and read the time on
its first line:

```
sbatch --array=1 scripts/slurm/scen1.slurm dec2024
grep "rep 1" slurm_logs/scen1_<jobid>_1.out       # rep 1 done in <N> s
```

Set `BLOCK` in that `.slurm` file to about `3600 / N`. Task 1 can keep running,
because the full run skips replicates that already have output.

**3. Full run**, with `--array=1-ceil(1000/BLOCK)`:

| Run | Command |
|---|---|
| scen0 | `sbatch --array=1-17 scripts/slurm/scen0.slurm dec2024` |
| scen1 | `sbatch --array=1-30 scripts/slurm/scen1.slurm oct2026` (BLOCK 34, sized on oct2026) |
| scen2 | `sbatch --array=1-30 scripts/slurm/scen2.slurm oct2026` (BLOCK 34, sized on oct2026) |
| scen3 | `sbatch --array=1-19 scripts/slurm/scen3.slurm dec2024` |
| scen3_strata | `sbatch --array=1-40 scripts/slurm/scen3_strata.slurm oct2026` |
| scen3_srs | `sbatch --array=1-34 scripts/slurm/scen3_srs.slurm oct2026` |
| scen3_unw | `sbatch --array=1-17 scripts/slurm/scen3_unw.slurm oct2026` |
| scen3_s2 | `sbatch --array=1-38 scripts/slurm/scen3_s2.slurm oct2026` |
| scen3_vunw | `sbatch --array=1-28 scripts/slurm/scen3_vunw.slurm none oct2026` |
| scen3_raw | `sbatch --array=1-50 scripts/slurm/scen3_raw.slurm oct2026` (10 CPUs per task, 10 reps in parallel) |
| scen4 | `sbatch --array=1-77 scripts/slurm/scen4_mi.slurm 1 oct2026` |
| scen5 | `sbatch --array=1-77 scripts/slurm/scen5_mi.slurm 1 oct2026` |

For the MI scenarios, the first argument is the variant (1–3). For scen3_vunw, it is the
step 2 weighting (`none` or `hh`), and each one is a separate array.

**4. Coverage**, once the estimates exist:

```
sbatch scripts/slurm/submit_coverage.slurm
```

This writes one table per run to `data/summaries/<pop>/coverage_<scen>.csv`, with
bias, empirical vs. model SE, model-based coverage, and coverage from the empirical SE.
