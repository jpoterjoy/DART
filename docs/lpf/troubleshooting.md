# LPF troubleshooting

[← LPF overview](../../README_LPF.md) · [Algorithm](algorithm.md) · [Namelist](namelist.md) · [Tuning](tuning.md)

## The run does not know the `pf_*` options

`filter` stops while reading `&assim_tools_nml`, or reports an unknown
namelist entry such as `frac_neff` or `pf_maxiter`. This means it was built
with the stock `assim_tools_mod.f90`. The LPF lives in
`assim_tools_mod.pf.f90`, which DART's build skips (it is listed as a
"nuisance" file in `build_templates/buildfunctions.sh`). To fix it:

```bash
cd assimilation_code/modules/assimilation
cp assim_tools_mod.pf.f90 assim_tools_mod.f90
cd ../../../models/<your_model>/work && ./quickbuild.sh
```

## Adaptive inflation is on

The LPF does not support adaptive prior or posterior inflation, which relies
on EAKF quantities. Set `inf_flavor = 0, 0` in `&filter_nml`. Fixed prior
inflation is applied outside the LPF and can be used, though it was at most
marginally helpful in tests. If the ensemble collapses, tune `frac_neff`, the
number of passes and `cutoff` first (see [Tuning](tuning.md)).

## `WARNING! Neff is ... but target is ...`

This comes from `pf_regularization`: for one variable, no tempering
coefficient $\beta$ reached the target effective ensemble size. The variable
then gets $\beta = 0$ and is not updated in that pass. Occasional warnings are
harmless. Frequent ones point to extremely sharp weights. Check for:
- observation error variances that are too small;
- gross outliers that pass QC (`outlier_threshold`);
- `frac_neff` set very close to 1.

## Every cycle uses `pf_maxiter` passes

`Number of PF iterations` equal to `pf_maxiter` in every cycle means the
tempering never finishes, so part of each observation's information is
dropped. Raise `pf_maxiter`, or lower `frac_neff` so each pass takes a bigger
step. With `pf_diag = .true.`, the `PFDIAG` lines show the mean $\beta$ and
the fraction of variables still tempered ($\beta \lt  1$) in each pass.

## The ensemble is over-dispersive (spread larger than error)

This is common with a single pass at small $N_e$. In bgrid at $N_e = 80$ the
spread was 2.1 against an RMSE of 1.35. It improves with tempering
(`pf_maxiter > 1`) and with ensemble size. Very high `frac_neff` (0.9 and
above) also inflates spread, because it forces many small passes.

## The filter diverges or produces NaN with a nonlinear observation operator

The LPF updates the observation-space priors by the same weighting and merging
as the state, and does **not** recompute $h(x)$ between passes. With a strongly
nonlinear $h$ and many passes, the observation priors can drift away from
$h$ of the merged state. An observation then stops constraining outlying
particles, which can run away in the forecast. In tests with the MATLAB
version of the LPF, this drift was the cause of such failures.

Mitigations:
- use fewer passes (lower `frac_neff`, or a smaller `pf_maxiter`);
- keep KDDM on (`pf_kddm = 1`);
- check that observations far outside the prior are being rejected by QC.

## Runs are slow

First set `pf_profile = .true.` for a few cycles: the `PFPROF` lines show how
much time goes to the passes, the exchanges between tasks and KDDM, and how
much to the rest of `filter`. Then, in order of payoff:
1. Compile with `-O3` and an architecture flag (e.g.
   `-march=skylake-avx512`) so that `exp` and `log` are vectorized. This
   roughly halves the run time.
2. Use binning (`pf_bin_rho = 0.95`): 15–27% faster in bgrid, and usually
   more accurate.
3. Keep `pf_cache_close = .true.`.
4. Reduce the number of passes (`frac_neff`, `pf_maxiter`).

The KDDM step is a small part of the cost.

## Memory

`pf_cache_close` stores about $N_e$ reals for every observation-variable pair
within the localization radius, on each task. With dense observations and a
wide `cutoff`, set it to `.false.` if memory is tight; runs are then a few
percent slower.

The exchanges of per-observation values between tasks use up to
[`pf_gather_mb`](namelist.md#pf_gather_mb) MB per task (default 200). With
very many observations per cycle and little memory per task, lower it; results
don't change, only the number of exchanges.

## Results change with the number of MPI tasks

The resampling random seed includes the MPI task number, so the same
experiment on a different number of tasks gives a different (equally valid)
realization. Results are reproducible for a fixed task count and `pf_seed`.
To compare configurations, use the same task count and several values of
`pf_seed`.
