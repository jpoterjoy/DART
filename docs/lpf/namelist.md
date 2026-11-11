# LPF namelist reference

[← LPF overview](../../README_LPF.md) · [Algorithm](algorithm.md) · [Tuning](tuning.md) · [Troubleshooting](troubleshooting.md)

All options are in `&assim_tools_nml`. The LPF-specific ones are listed
first, then the standard DART options that interact with the LPF. Defaults
are those in `assim_tools_mod.pf.f90`.

**LPF options:** [`filter_kind`](#filter_kind) ·
[`cutoff`](#cutoff) · [`frac_neff`](#frac_neff) · [`pf_maxiter`](#pf_maxiter) ·
[`pf_kddm`](#pf_kddm) · [`bw_factor`](#bw_factor) · [`pf_merge`](#pf_merge) ·
[`pf_nmoments`](#pf_nmoments) · [`pf_qm_moments`](#pf_qm_moments) ·
[`pf_bin_rho`](#pf_bin_rho) · [`pf_cache_close`](#pf_cache_close) ·
[`pf_seed`](#pf_seed) · [`pf_diag`](#pf_diag) ·
[`pf_gather_mb`](#pf_gather_mb) · [`pf_profile`](#pf_profile) ·
[`pf_enkf_hybrid`](#pf_enkf_hybrid) · [`min_residual`](#min_residual) ·
[`pf_kf_rtps_coeff`](#pf_kf_rtps_coeff) · [`pf_alpha`](#pf_alpha)

**Also:** [vertical localization tables](#vertical-localization-tables-cov_cutoff_nml)
for observations without a vertical coordinate (`&cov_cutoff_nml`).

---

### `filter_kind`
**integer, default `9`**

`9` selects the LPF. The other values are DART's standard filters (1 EAKF,
2 EnKF, 3 kernel, 4 particle, 5 random draw, 6 deterministic draw, 8 RHF);
none of the LPF options below affect them.

### `cutoff`
**real, default `0.3`**

Localization half-width (Gaspari–Cohn), in the location module's units
(radians for `threed_sphere`). An observation affects variables out to
`2*cutoff`. In the LPF it sets $\rho$ for both the localized weights and the
merge. Per-type values: `special_localization_obs_types` and
`special_localization_cutoffs`.

### `frac_neff`
**real, default `0.80`**

Target effective ensemble size for tempering, as a fraction of $N_e$. Each
variable's $\beta$ is the largest that keeps
$N_\text{eff} \ge$ `frac_neff` × $N_e$.
- **Higher** means smaller tempering steps, and more passes when
  `pf_maxiter > 1`.
- **Lower** means bigger steps and fewer passes, at greater risk of weight
  collapse.

Range 0 < `frac_neff` < 1. See the [tuning guide](tuning.md#frac_neff-and-the-number-of-passes).

### `pf_maxiter`
**integer, default `1`**

Maximum number of LPF passes per cycle.
- `1`: a single pass; any tempering residual left after it is discarded.
  The same happens to any residual left when a tempered run reaches
  `pf_maxiter`.
- `> 1`: tempered LPF, which stops early once every observation's tempering
  is used up.

The log prints `Number of PF iterations` each cycle. A run that always hits
`pf_maxiter` is not finishing its tempering.

### `pf_kddm`
**integer, default `0`**

Kernel density distribution mapping after each pass, for variables whose
particles moved outside the prior range:
- `0`: off.
- `1`: KDDM, extrapolating linearly beyond the particle range.
- `2`: bounded KDDM, which keeps particles within the range of the prior
  particles.

### `bw_factor`
**real, default `1.0`**

Multiplier on the KDDM kernel bandwidth $\sqrt2~\sigma~N_e^{-1/5}$. Used
only with `pf_kddm > 0`.

### `pf_merge`
**integer, default `0`**

How resampled and current particles are merged:
- `0`: moment matching with a Cornish-Fisher transform (`pf_update`;
  moments set by `pf_nmoments`).
- `1`: rank quantile mapping onto the weighted posterior (`pf_update_qm`;
  adjustment set by `pf_qm_moments`).

### `pf_nmoments`
**integer, default `4`**

Moments fitted by the moment-matching merge (`pf_merge = 0`):
- `2`: mean and variance.
- `3`: adds skewness.
- `4`: adds kurtosis.

Higher moments are shrunk by their sampling reliability, so `4` is safe at
small $N_e$.

### `pf_qm_moments`
**integer, default `2`**

With `pf_merge = 1`, the adjustment after quantile mapping:
- `2`: shift and scale to the exact posterior mean and variance.
- `1`: shift to the mean only.
- `0`: none.

### `pf_bin_rho`
**real, default `0.0` (off)**

Merge-only binning of nearby observations. An observation joins a bin when
its localization coefficient with the bin centroid is at least `pf_bin_rho`.
Resampling and merging then happen once per bin, while the weights still use
each observation's own localization.
- `0`: off.
- `1`: bin only co-located observations.
- `0.95` is a good starting value (members within about 0.18 × cutoff of the
  centroid).

### `pf_cache_close`
**logical, default `.true.`**

Keep each observation's close-variable lists and localization factors across
passes, and its log-weight increments between the regularization and update
passes. This makes passes cheaper. It costs memory of about $N_e$ reals per
observation-variable pair on each task. It is turned off automatically with
adaptive localization or localization diagnostics.

### `pf_seed`
**integer, default `0`**

Offset for the resampling random seed. Different values give independent
realizations of the same experiment.

### `pf_diag`
**logical, default `.false.`**

Print one `PFDIAG` line per pass, covering state variables with observations
nearby:
- the mean $\beta$;
- the fraction with $\beta \lt  1$;
- the mean untempered $N_\text{eff}/N_e$.

### `pf_gather_mb`
**real, default `200.0`**

Memory per MPI task, in MB, for exchanging per-observation values between
tasks: weights, QC and verticals, and resampling indices. These are exchanged
in blocks of observations, one exchange per block, rather than one broadcast
per observation. Each observation takes about $3N_e + 4$ reals, so this sets
the number of observations per block. When all observations fit in one block,
each pass exchanges them once.

Results do not depend on it. Lower it only if memory is tight with very many
observations per cycle.

### `pf_profile`
**logical, default `.false.`**

Print wall-clock timings of the LPF, on task 0:
- `PFPROF iter` and `PFPROF2 iter` lines for each pass: setup, the
  regularization loop and its exchanges, the tempering coefficients, the
  sequential loop and its exchanges, KDDM, and, within the sequential loop,
  the owner's work, the state and obs-space updates and the bin merges. Each
  is given for task 0 and as the maximum over tasks.
- A `PFPROF call` line per cycle: the time in the LPF passes, and the time
  since the previous cycle's LPF ended (model advance, I/O, forward
  operators).

Use it to find where time goes before tuning for speed; it doesn't change
results.

### `pf_enkf_hybrid`
**logical, default `.false.`**

Finish each cycle with an EAKF pass (hybrid PF-EAKF, Poterjoy 2021). The PF
passes then assimilate only 1 − `min_residual` of each
observation's information.

### `min_residual`
**real, default `0.0`**

Hybrid only. The fraction of each observation's information left for the
EAKF pass, which uses error variance $R$ / `min_residual`. It is
ignored when `pf_enkf_hybrid = .false.`.

### `pf_kf_rtps_coeff`
**real, default `0.0`**

Hybrid only. Relaxation-to-prior-spread coefficient applied after the EAKF
pass (0 = none, 1 = full prior spread).

### `pf_alpha`
**real, default `1.0`: deprecated**

Not used. It is kept so that existing namelists still read.

---

## Standard DART options that interact with the LPF

| Option | Default | With the LPF |
|---|---|---|
| `special_localization_obs_types`, `special_localization_cutoffs` | none | Per-type cutoffs for weights, merges and bins. |
| `adjust_obs_impact`, `obs_impact_filename`, `allow_any_impact_values` | off | Impact factors multiply $\rho$ in both passes. |
| `adaptive_localization_threshold`, `adaptive_cutoff_floor` | `-1`, `0` | Supported; turns off `pf_cache_close`. |
| `output_localization_diagnostics`, `localization_diagnostics_file` | off | Supported; turns off `pf_cache_close`. |
| `convert_all_obs_verticals_first`, `convert_all_state_verticals_first` | `.true.`, `.false.` | Both passes use the converted observation vertical. |
| `distribute_mean`, `close_obs_caching`, `print_every_nth_obs` | `.false.`, `.true.`, `0` | As in DART. |
| `sort_obs_inc`, `spread_restoration`, `sampling_error_correction`, `rectangular_quadrature`, `gaussian_likelihood_tails` | off | Not used by the LPF. |

In `&filter_nml`, do not use adaptive inflation with the LPF. It relies on
EAKF quantities the LPF does not compute. Fixed (non-adaptive) prior inflation
is applied by `filter` outside the LPF and does work: in one Lorenz-96 test a
fixed factor of 1.05 helped marginally. Usually leave it off
(`inf_flavor = 0, 0`).

---

## Vertical localization tables (`&cov_cutoff_nml`)

Some observations have no single vertical location: a satellite channel, or a
layer-mean temperature, depends on several model levels at once. Located with
`VERTISUNDEF`, such an ob gets horizontal localization only, so it updates
every level fully. These tables give it a vertical localization factor
instead, read from the namelist for each obs type. They live in DART's
`cov_cutoff_mod` and apply to every DART filter, not only the LPF, and to any
model.

| Option | Type, default | Meaning |
|---|---|---|
| `vloc_types` | strings, none | Obs types that use a table (up to 20). |
| `vloc_nlev` | integer, `0` | Number of rows in every table (1–200). Needed when `vloc_types` is set. |
| `vloc_factors(k, t)` | real, `1.0` | Factor in row `k` of the table for the `t`-th type in `vloc_types`. |
| `vloc_coord(t)` | string, `'level'` | The vertical coordinate of table `t`: `'level'`, `'pressure'`, `'height'` or `'scale_height'`. |
| `vloc_points(k, t)` | real, unset | For non-level tables: the coordinate of row `k` (Pa, m, or scale height). Must be set for rows 1 to `vloc_nlev`, strictly increasing or decreasing; positive for pressure. Not used by level tables. |
| `vloc_surface_level` | integer, `1` | Row whose factor applies to surface (`VERTISSURFACE`) targets, in any table. |
| `vloc_mismatch` | string, `'error'` | What happens when a table meets a target located in a different vertical coordinate: `'error'` stops the run; `'warn'` uses a factor of 1, with one warning per table and coordinate. |

How a factor is applied, for an ob of a listed type located with `VERTISUNDEF`:
- **Level tables:** row `k` is model level `k`. The target's level is
  interpolated linearly between rows.
- **Pressure, height and scale-height tables:** the target's vertical
  coordinate is interpolated between the rows' `vloc_points`, linearly in
  log-pressure for pressure tables, linearly otherwise. Beyond the first or
  last row, that row's factor is used.
- **Surface targets** use row `vloc_surface_level`; **targets without a
  vertical** get 1.
- **A target in any other coordinate** is handled by `vloc_mismatch`. Before
  this check was added, such targets silently got a factor of 1.
- Obs types not listed, and listed types located *with* a vertical, are
  unchanged.
- In the LPF, an ob with a table is binned only with obs of the same type.

**Match the table to how your model locates its state.** The table's
coordinate must be the one the *targets* use when the localization is
computed: the state variables, and any other obs the ob updates. bgrid
locates its state on model levels, so level tables work there. Many models
(e.g. WRF, MPAS, CAM) localize in pressure or height, often after converting
the state verticals (`convert_all_state_verticals_first`, or the model's own
vertical-localization setting). For those, use a `'pressure'` or `'height'`
table. Pressure tables also suit channel weighting functions, which are
usually given in pressure. The default `vloc_mismatch = 'error'` stops the run
at the first target in the wrong coordinate, so a table can't be ignored
unnoticed.

**Which factors to use: the support rule.** Set the factor to 1 on every level
the ob depends on (where its weighting function is nonzero), and taper with
Gaspari–Cohn outside that range. In bgrid tests with three synthetic
layer-temperature channels (6 seeds, prior RMSE, days 60–120):

| Localization of the channels | u | v | T | ps |
|---|---|---|---|---|
| **Support rule** | **2.174** | **1.913** | **1.217** | 191.7 |
| Factor 1 on the main lobe, Gaspari–Cohn outside | 2.186 | 1.925 | 1.562 | 190.2 |
| Weighting function smoothed with Gaspari–Cohn | 2.180 | 1.918 | 1.583 | 189.4 |
| Horizontal only (no table) | 2.276 | 2.029 | 1.701 | 195.5 |
| EAKF (1 seed) | 2.244 | 1.991 | 1.273 | 210.9 |

The support rule cut T error by 22% compared with the main-lobe rule, and beat
the EAKF on all four variables.
Two schemes that damp the update on levels where the channel *does* depend
(the factor at the weighting-function peak, and the weighting function
normalized by its maximum) made the tempered LPF diverge. On a network of
channels with long tails it tied with the main-lobe rule. On eight narrow
channels it was 1–1.5% worse than the main-lobe rule.

Examples (the type names are illustrative).

A level table, for a 5-level model and a channel that depends on levels 2–4:

```
&cov_cutoff_nml
   vloc_types          = 'LAYER_TEMPERATURE_CH1',
   vloc_nlev           = 5,
   vloc_surface_level  = 5,
   vloc_factors(1:5,1) = 0.685, 1.0, 1.0, 1.0, 0.685,
/
```

A pressure table, for a channel that depends on 850–300 hPa in a model that
localizes in pressure. Rows run from the surface up; the factor is 1 across
the support and tapers above and below it:

```
&cov_cutoff_nml
   vloc_types          = 'SATELLITE_CHANNEL_7',
   vloc_coord          = 'pressure',
   vloc_nlev           = 6,
   vloc_surface_level  = 1,
   vloc_points(1:6,1)  = 100000.0, 85000.0, 60000.0, 30000.0, 20000.0, 10000.0,
   vloc_factors(1:6,1) =      0.6,     1.0,     1.0,     1.0,     0.6,     0.1,
/
```
