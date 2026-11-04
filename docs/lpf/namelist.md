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
[`pf_enkf_hybrid`](#pf_enkf_hybrid) · [`min_residual`](#min_residual) ·
[`pf_kf_rtps_coeff`](#pf_kf_rtps_coeff) · [`pf_alpha`](#pf_alpha)

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
$N_\text{eff} \ge \texttt{frac\_neff}\cdot N_e$.
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

Multiplier on the KDDM kernel bandwidth $\sqrt2\,\sigma\,N_e^{-1/5}$. Used
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
- the fraction with $\beta < 1$;
- the mean untempered $N_\text{eff}/N_e$.

### `pf_enkf_hybrid`
**logical, default `.false.`**

Finish each cycle with an EAKF pass (hybrid PF-EAKF, Poterjoy 2021). The PF
passes then assimilate only $1 - \texttt{min\_residual}$ of each
observation's information.

### `min_residual`
**real, default `0.0`**

Hybrid only. The fraction of each observation's information left for the
EAKF pass, which uses error variance $R/\texttt{min\_residual}$. It is
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
