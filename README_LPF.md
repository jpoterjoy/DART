# Localized particle filter (LPF) in DART

This branch carries the localized particle filter (LPF) of Poterjoy (2016),
Poterjoy and Anderson (2016), Poterjoy et al. (2019) and Poterjoy (2021),
selected with `filter_kind = 9` in `&assim_tools_nml`. This document describes
how the Fortran implementation works and what each LPF namelist option does.

Questions: Jon Poterjoy (poterjoy@umd.edu).

## Where the code is

| File | Role |
|---|---|
| `assimilation_code/modules/assimilation/assim_tools_mod.pf.f90` | The LPF version of `assim_tools_mod` (all LPF code lives here). |
| `assimilation_code/modules/assimilation/assim_tools_mod.f90` | The module DART compiles. Stock DART in git. |

`assim_tools_mod.pf.f90` is listed as a "nuisance" file in
`build_templates/buildfunctions.sh`, so it is **not** compiled as-is. To build
filter with the LPF, copy it over the compiled module before building:

```
cd assimilation_code/modules/assimilation
cp assim_tools_mod.pf.f90 assim_tools_mod.f90
cd ../../../models/<model>/work
./quickbuild.sh            # or ./quickbuild.sh nompi filter
```

Restore the stock module with `git checkout assim_tools_mod.f90` when done.

The LPF supports only `filter_kind = 9` (with an optional EAKF final step, see
the hybrid below). Do not combine it with adaptive prior or posterior
inflation.

## How the LPF works

The ensemble of `Ne` members is treated as particles. Observations are
assimilated one at a time, as in DART's other filters, but each observation
updates the particles by weighting, resampling and merging rather than by a
Kalman-type increment. Each state variable and each observation prior carries
its own set of particle weights, which accumulate the localized influence of
all observations processed so far.

### Localized weights

For observation `y` with error variance `R` and prior obs-space particles
`h_i`, the (untempered) likelihood weights are

```
w_i  ∝ exp( -(y - h_i)^2 / (2R) ),     sum_i w_i = 1.
```

For a state variable (or another observation's prior) at localization
coefficient `ρ` from the observation (Gaspari-Cohn of the distance, with
`cutoff` as the half-width), the localized weight is

```
ω_i = 1 + ρ (Ne w_i - 1),
```

which equals the full weight at `ρ = 1` and is uniform at `ρ = 0`. Each
variable keeps the running sum of `-log ω_i` over the observations
processed so far (`lw` for state variables, `lhw` for observation priors).

**Weight floor.** Before forming `w`, the code checks whether the smallest
normalized weight falls below `1e-15/Ne`. If it does, it tempers the
likelihood (an effective inflation of `R`, found by bisection in
`pf_regularization_minw`) until the smallest weight reaches that floor. This
is hard-coded and acts as a safeguard for very sharp likelihoods.

### Tempering (regularization) and iterations

Assimilating all observations at once can collapse the weights. The LPF
therefore tempers each variable's accumulated weights with a coefficient
`β ∈ (0, 1]`: the weights `exp(-β · lw)` are required to have an effective
ensemble size of at least `frac_neff · Ne`. Each variable gets its own `β`,
found in `pf_regularization` (safeguarded Newton iteration on
`Neff(β) = frac_neff · Ne`; `β = 1` if the untempered weights already meet
the target).

One assimilation cycle can then take several passes (iterations) over all
observations:

1. **Regularization pass.** Compute every observation's weights and every
   variable's accumulated `-log` localized weights, then each variable's `β`
   (state space) and each observation's `β` (obs space). Each variable keeps a
   residual `1 - Σ β` of tempering not yet used.
2. **Update pass.** Assimilate the observations sequentially with the tempered
   weights (see below), then optionally apply KDDM (see `pf_kddm`).
3. Repeat with the updated particles until every observation's residual is
   used up (`Σ β = 1`) or `pf_maxiter` passes have been made.

With `pf_maxiter = 1` the filter is a single pass (`β` still limits weight
collapse, but any remaining residual is dropped). Larger values give a
tempered, multi-pass LPF. The log reports `Number of PF iterations` each
cycle.

### Update pass (one observation)

For each observation in sequence:

1. **Resampling.** The owner of the observation draws `Ne` indices by
   systematic resampling (`pf_sample`) from the observation's own tempered,
   localized weights. Surviving particles keep their indices (to avoid
   needless shocks), and duplicates fill the vacated slots in random order.
   If the weights are nearly uniform, resampling is skipped.
2. **Weights for each nearby variable.** Each state variable and observation
   prior within the localization radius adds this observation's `-log ω` to
   its running sum, and forms its posterior weights `exp(-β · lw)`. It then
   computes the weighted posterior mean and variance from its prior
   particles.
3. **Merge.** The resampled particles (prior values at the drawn indices) are
   merged with the current particles. Their perturbations about the posterior
   mean are blended with weight `β·ρ`:

   ```
   x_blend = m + βρ (x_resampled - m) + (1 - βρ) (x_current - m).
   ```

   The blended particles are then adjusted to match the weighted posterior.
   Two methods are available (see `pf_merge`):
   - **Moment matching** (default, `pf_update`). Fit up to four moments (mean,
     variance, skewness, kurtosis; `pf_nmoments`) with a Cornish-Fisher
     transform. The skewness and kurtosis corrections are shrunk by their
     sampling reliability and clipped to keep the map near monotone, then the
     result is rescaled to the exact posterior mean and variance.
   - **Rank quantile mapping** (`pf_update_qm`). Keep only the ranks of the
     blend and assign each particle the weighted posterior quantile at
     `(k - 1/2)/Ne`, with short kernel tails beyond the particle range.
4. Observation priors are updated the same way, including observations
   already assimilated, because later passes reuse them.

Note: resampling happens at the observation location, but the merge at each
variable is localized (`βρ`). A variable far from the observation keeps its
own particles, and a variable at the observation is fully replaced by
resampled particles.

### Optional steps

- **KDDM** (`pf_kddm > 0`). After each pass, every state variable and
  observation prior whose posterior particles fall outside the span of its
  prior particles is mapped by kernel density distribution mapping
  (`pf_kddm_update`). Prior particles are mapped to the quantiles of the
  kernel-smoothed weighted posterior, which removes spurious spread outside
  the prior range.
- **Binning** (`pf_bin_rho > 0`). Nearby observations are grouped into bins
  (`pf_make_bins`). Within a bin, each observation still contributes its own
  localized weights, but resampling and merging happen once, after the bin's
  last observation, localized about the bin centroid. This avoids repeated
  resampling for clustered observations and is faster.
- **Hybrid PF-EAKF** (`pf_enkf_hybrid = .true.`). The PF passes stop when
  each observation has used `1 - min_residual` of its tempering. A final pass
  then assimilates the observations with the EAKF, with the observation
  error variance inflated by `1/min_residual`. Relaxation to prior spread
  (`pf_kf_rtps_coeff`) is applied after the EAKF step.

### Random numbers

The random sequence used in resampling is initialized once per run. Its seed
depends on the time of the first observation, the MPI task and `pf_seed`, so
different cycles and iterations do not repeat the same random numbers. Vary
`pf_seed` to get independent realizations of the same experiment.

## Namelist options (`&assim_tools_nml`)

### LPF options

| Option | Type, default | Meaning |
|---|---|---|
| `filter_kind` | integer, `9` | `9` selects the LPF. (Other values are DART's standard filters: 1 EAKF, 2 EnKF, 3 kernel, 4 particle, 5 random draw, 6 deterministic draw, 8 RHF; the LPF code is not used for those.) |
| `cutoff` | real, `0.3` | Localization half-width (Gaspari-Cohn), in the location module's units (radians for `threed_sphere`). Sets `ρ` for the localized weights and the merge. Per-type values via `special_localization_obs_types` / `special_localization_cutoffs`. |
| `frac_neff` | real, `0.80` | Target effective ensemble size for tempering, as a fraction of `Ne`. Higher values mean smaller steps, more tempering and (with `pf_maxiter > 1`) more passes. |
| `pf_maxiter` | integer, `1` | Maximum number of LPF passes per cycle. `1`: single pass. `> 1`: tempered, multi-pass LPF; passes stop early when all tempering is used up. |
| `pf_kddm` | integer, `0` | KDDM after each pass: `0` off, `1` KDDM with extrapolation beyond the particle range, `2` bounded KDDM (quantiles clipped to the posterior particle range). |
| `bw_factor` | real, `1.0` | Multiplier on the KDDM kernel bandwidth (`√2 · σ · Ne^(-1/5)`). |
| `pf_nmoments` | integer, `4` | Moments fitted by the default merge (`pf_merge = 0`): `2` mean and variance, `3` adds skewness, `4` adds kurtosis. |
| `pf_merge` | integer, `0` | Merge method: `0` moment matching (Cornish-Fisher, `pf_update`), `1` rank quantile mapping (`pf_update_qm`). |
| `pf_qm_moments` | integer, `2` | With `pf_merge = 1`: after the quantile mapping, `2` shift and scale to the exact posterior mean and variance, `1` shift to the mean only, `0` leave as mapped. |
| `pf_bin_rho` | real, `0.0` | Merge-only binning of nearby observations. An observation joins a bin when its localization coefficient with the bin centroid is at least `pf_bin_rho`. `0` disables binning; `1` bins only co-located observations. Typical value: `0.95`. |
| `pf_cache_close` | logical, `.true.` | Keep each observation's close-variable lists and localization factors across passes, and its log-weight increments between the regularization and update passes. Faster; costs about `Ne` reals per observation-variable pair on each task. Ignored with adaptive localization or localization diagnostics. |
| `pf_seed` | integer, `0` | Offset for the resampling random seed. Change it to get independent runs. |
| `pf_diag` | logical, `.false.` | Print, for each pass, the mean `β` over state variables with observations nearby, the fraction with `β < 1`, and their mean untempered effective ensemble size (`PFDIAG` lines in the log). |
| `pf_enkf_hybrid` | logical, `.false.` | Hybrid PF-EAKF: finish each cycle with an EAKF pass (see above). |
| `min_residual` | real, `0.0` | Hybrid only: the fraction of each observation's information left to the EAKF step (its error variance is inflated by `1/min_residual`). Ignored without the hybrid. |
| `pf_kf_rtps_coeff` | real, `0.0` | Hybrid only: relaxation-to-prior-spread coefficient applied after the EAKF step. |
| `pf_alpha` | real, `1.0` | **Deprecated, not used.** Kept so existing namelists still read. |

### Standard DART options relevant to the LPF

These options are inherited from DART. The ones that interact with the LPF:

| Option | Default | Note for the LPF |
|---|---|---|
| `special_localization_obs_types`, `special_localization_cutoffs` | none | Per-type cutoffs, used for weights, merges and bins. |
| `adjust_obs_impact`, `obs_impact_filename`, `allow_any_impact_values` | off | Impact factors multiply `ρ` in both passes. |
| `adaptive_localization_threshold`, `adaptive_cutoff_floor` | `-1`, `0` | Adaptive localization works, but turns off `pf_cache_close`. |
| `output_localization_diagnostics`, `localization_diagnostics_file` | off | Turns off `pf_cache_close`. |
| `convert_all_obs_verticals_first`, `convert_all_state_verticals_first` | `.true.`, `.false.` | Vertical conversion for localization; both passes use the converted vertical. |
| `distribute_mean`, `close_obs_caching` | `.false.`, `.true.` | As in DART. |
| `print_every_nth_obs` | `0` | Progress messages. |
| `sort_obs_inc`, `spread_restoration`, `sampling_error_correction`, `rectangular_quadrature`, `gaussian_likelihood_tails` | off | Not used by the LPF (EAKF/RHF options). |

Adaptive inflation (`&filter_nml` `inf_flavor`) should be off with the LPF.

## Example

A tempered LPF with up to 3 passes, KDDM, four-moment merge and binning:

```
&assim_tools_nml
   filter_kind      = 9,
   cutoff           = 0.5,
   frac_neff        = 0.8,
   pf_maxiter       = 3,
   pf_kddm          = 1,
   pf_nmoments      = 4,
   pf_merge         = 0,
   pf_bin_rho       = 0.95,
   pf_cache_close   = .true.,
   pf_seed          = 0,
   /
```

A single-pass LPF: the same with `pf_maxiter = 1`.

## Practical notes

- **Tuning.** `frac_neff` and `cutoff` matter most. With few, short-interval
  observations, `frac_neff ≈ 0.8` with 2–3 passes works well. Long intervals
  or sharp likelihoods need lower `frac_neff` and more passes. In tests, the
  number of passes mattered more than `frac_neff` itself.
- **Cost.** Each pass costs about one EAKF pass plus the regularization pass.
  Most of the time goes into `exp`/`log` in the weight loops. Build with
  `-O3` and an architecture flag (e.g. `-march=skylake-avx512`) to use
  vectorized math functions; this roughly halves the run time.
- **Binning** (`pf_bin_rho = 0.95`) helped and was faster in every Lorenz-96
  and bgrid configuration tested.
- **The weight floor** (`1e-15/Ne`) is deliberate. Removing it degraded
  results for sharp likelihoods.

## References

- Poterjoy, J., 2016: A localized particle filter for high-dimensional
  nonlinear systems. *Mon. Wea. Rev.*, 144, 59–76.
- Poterjoy, J., and J. L. Anderson, 2016: Efficient assimilation of simulated
  observations in a high-dimensional geophysical system using a localized
  particle filter. *Mon. Wea. Rev.*, 144, 2007–2020.
- Poterjoy, J., L. Wicker, and M. Buehner, 2019: Progress toward the
  application of a localized particle filter for numerical weather
  prediction. *Mon. Wea. Rev.*, 147, 1107–1126.
- Poterjoy, J., 2021: Implications of multivariate non-Gaussian data
  assimilation for multiscale weather prediction. *Mon. Wea. Rev.* (hybrid
  PF-EAKF).
