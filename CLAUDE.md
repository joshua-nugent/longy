# longy — Claude Code Context

## What is longy?

longy ("lawn-gee") is a private R package for longitudinal causal inference — an
alternative to `ltmle` and `stremr`. It provides IPW, G-computation, and TMLE
estimators with IC/EIF/bootstrap/sandwich inference, plus cross-fitting (CV-TMLE).

## Architecture

- **S3 classes**: `longy_data` (pipeline object), `longy_result` (estimation result), `longy_contrast` (treatment effect contrast), `longy_influence_diag` (influence/weight diagnostics)
- **Plain functions** internally (no R6)
- **data.table** internally for performance; accept data.frame at boundary
- **SuperLearner** for ML with glm fallback (SuperLearner is in Suggests, not Imports)
- **Long format only**: one row per person-time
- **Unified pipeline**: all `fit_*` and `estimate_*` accept/return `longy_data`; results accumulate in `obj$results`

## Two-level API

1. `longy()` — high-level wrapper, does everything in one call
2. Modular pipeline:
   ```r
   longy_data() |> define_regime() |> fit_treatment() |> fit_censoring() |>
     fit_observation() |> compute_weights() |> estimate_ipw()
   ```

### `longy()` key parameters

- `estimator`: `"ipw"` (default), `"gcomp"`, `"tmle"`, `"both"` (ipw+gcomp), `"all"` (ipw+gcomp+tmle), or any `c(...)` combination. Unadjusted estimates are always computed. **IPW is not supported with competing risks** — `longy()` auto-drops it with a warning; `estimate_ipw()` errors.
- `risk_set_treatment`: `"all"` (default) or `"followers"` — who trains treatment models
- `risk_set_outcome`: `"all"` (default) or `"followers"` — who trains outcome models (G-comp/TMLE)
- `outcome_range`: Numeric(2) for continuous outcome scaling in TMLE (NULL = empirical)
- `g_bounds`: Bounds for unstabilized cumulative g in TMLE clever covariate (default `c(0.01, 1)`). IPW uses `truncation`/`truncation_quantile` instead.
- `competing`: Column for binary competing event indicator (survival outcomes)
- `cross_fit` / `cross_fit_seed`: CV fold count and seed
- `use_ffSL`: Use future-factorial SuperLearner (parallel CV)
- `parallel`: Dispatch time-point models in parallel via `future.apply`
- `k`: Lag depth for covariate history (default `0`)
- `min_obs`, `min_events`, `adaptive_cv`: Model fitting controls

## DAG ordering assumption

longy assumes the within-period temporal ordering **L(t) → A(t) → C(t) → Y(t)**:

1. **L(t)**: Time-varying covariates are realized
2. **A(t)**: Treatment is assigned (may depend on L(t) and history)
3. **C(t)**: Censoring occurs (may depend on A(t), L(t), and history)
4. **Y(t)**: Outcome is measured (if uncensored and observed)

This means A(t) is a valid covariate for censoring and observation models, and is
included by default. This matches ltmle's vignette ordering (W → A1 → C → L → A2 → Y)
and lmtp's structure. Note that stremr assumes L-C-A-Y by default.

### Default covariates by model

| Model | Default covariates | Rationale |
|-------|-------------------|-----------|
| g_A (treatment) | baseline + timevarying | A(t) cannot condition on itself |
| g_C (censoring) | baseline + timevarying + **A(t)** | A(t) precedes C(t) in DAG |
| g_R (observation) | baseline + timevarying + **A(t)** | A(t) precedes R(t) in DAG |
| Q (outcome) | baseline + timevarying + **A(t)** | A(t) precedes Y(t) in DAG |

All models also include lag columns (`.longy_lag_*`) when `k > 0`.

## Three nuisance models

- **g_A** (treatment): P(A(t)=1 | past) — who gets treated?
- **g_C** (censoring): P(C(t)=0 | past) — who stays? (absorbing, cumulated)
- **g_R** (observation): P(R(t)=1 | past) — whose outcome is measured? (intermittent, NOT cumulated)

## Censoring interface

The `censoring` parameter accepts a **single character/factor column name**. The
column must contain `"uncensored"` for non-censored rows; any other value is
treated as a censoring cause.

- Single cause: `C = c("uncensored", "censored", ...)` — one internal indicator `.cens_censored`
- Multiple causes: `C = c("uncensored", "death", "ltfu", ...)` — creates `.cens_death` and `.cens_ltfu`

Internally, `longy_data()` decomposes the character column into binary `.cens_<cause>`
columns. These are stored in `nodes$censoring` (character vector of internal column names).
The original column name is in `nodes$censoring_col` and the cause labels in
`nodes$censoring_levels`. All downstream code (`fit_censoring`, `compute_weights`,
TMLE, cross-fitting) operates on the binary column names and is unchanged.

## Key data structure: longy_data object

```r
list(
  data       = <data.table keyed on (id, time)>,
  nodes      = list(id, time, outcome, treatment,
                    censoring,         # char vector of internal binary col names (e.g. ".cens_censored")
                    censoring_col,     # original column name (e.g. "C"), or NULL
                    censoring_levels,  # non-"uncensored" levels (e.g. c("death","ltfu")), or NULL
                    observation, cluster, sampling_weights,
                    baseline, timevarying, outcome_type, competing,
                    lag_vars,          # char vector of columns with lag columns created
                    lag_k),            # integer or Inf — lag depth
  regimes    = list(),
  fits       = list(treatment = list(), censoring = list(), observation = list(), outcome = list()),
  weights    = list(),
  results    = list(),   # {regime}_{estimator} -> longy_result (e.g. always_ipw, always_tmle)
  crossfit   = list(enabled = FALSE, n_folds = NULL, fold_id = NULL),
  meta       = list(n_subjects, n_obs, n_times, time_values, max_time, min_time)
)
```

### Refit protection

- `fit_*(refit=FALSE)`: errors if already fitted for requested regime. Use `refit=TRUE` to override.
- `compute_weights(recompute=FALSE)`: same pattern.
- `estimate_*`: no protection — silently overwrites (re-estimation is cheap).

### Result keys

Results stored as `obj$results$<regime>_<estimator>`:
- `always_ipw`, `always_gcomp`, `always_tmle`, `never_ipw`, etc.
- Access via `results(obj, regime="always", estimator="tmle")`

## File map

| File | Purpose |
|------|---------|
| R/data_input.R | `longy_data()` constructor, validation, censoring decomposition, `set_crossfit()` |
| R/regimes.R | `define_regime()` + `.evaluate_regime()` |
| R/fit_treatment.R | g_A models + `.add_tracking_columns()` / `.remove_tracking_columns()` |
| R/fit_censoring.R | g_C models (loops over `nodes$censoring` binary columns) |
| R/fit_observation.R | g_R models (intermittent) |
| R/fit_outcome.R | Outcome models via sequential regression (G-comp/TMLE) |
| R/weights.R | `compute_weights()` + `.compute_cumulative_g()` shared helper |
| R/estimate_ipw.R | `estimate_ipw()` — IPW estimator (not supported with competing risks) |
| R/contrast.R | `contrast()` — treatment effect contrasts between regimes + `longy_contrast` S3 class |
| R/longy_result.R | `print/summary/plot.longy_result` + `.estimator_label()` for all estimators |
| R/estimate_gcomp.R | `estimate_gcomp()` — G-comp estimator + bootstrap |
| R/estimate_tmle.R | `estimate_tmle()` — TMLE estimator + EIF inference |
| R/estimate_unadjusted.R | `estimate_unadjusted()` — naive means among regime-followers (reference) |
| R/learner-adaptations.R | Docs + wrapper functions for SL learner compatibility (SL.xgboost.reg, SL.glmnet.reg, SL.ranger.longy) |
| R/inference.R | IC-based + bootstrap inference (IPW, G-comp, TMLE) |
| R/diagnostics.R | `weight_diagnostics()`, `positivity_diagnostics()`, `prediction_diagnostics()`, `sl_diagnostics()`, `clip_diagnostics()`, `tmle_diagnostics()`, `influence_diagnostics()`, `plot_sl_diagnostics()`, `plot_influence_diagnostics()`, `longy_influence_diag` S3 class |
| R/crossfit.R | Cross-fitting: `.cf_fit_treatment/censoring/observation()`, `.cf_estimate_tmle()` |
| R/longy.R | High-level `longy()` wrapper + `results()` accessor + `plot.longy_data` |
| R/ffSL.R | Future-factorial SuperLearner (parallel CV via `future.apply`) |
| R/utils.R | Internal helpers: `.as_longy_data()`, `.safe_sl()`, `.bound()`, `.predict_from_fit()`, `.vmsg()`, `.get_lag_covariates()` |
| R/longy-package.R | Package-level docs, `@import data.table`, `globalVariables()`, `sim_longy` docs |
| R/zzz.R | `.onLoad()` hook |

## Coding conventions

- snake_case for everything
- Internal helpers prefixed with `.` (e.g., `.bound()`, `.safe_sl()`)
- data.table internally, data.frame at boundaries
- Use `list()` instead of `.()` in data.table `[` calls (avoids cedta() errors)
- All exported functions get roxygen2 docs
- Tests use testthat edition 3
- SuperLearner access via `requireNamespace("SuperLearner", quietly = TRUE)`

## Weight computation (critical logic)

Predicted probabilities are stored unbounded by `fit_treatment`, `fit_censoring`,
and `fit_observation`. Bounding happens in `compute_weights()`, which bounds the
**product** of raw probabilities (g_a * g_c * g_r) before computing weights.
The bounding adjustment is absorbed into the AC (absorbing) component so that
observation (intermittent) remains point-in-time.

```
g_a(t)      = P(A=d(t)|past)           # from fit_treatment (raw, unbounded)
g_c(t)      = P(C=0|past)              # from fit_censoring (raw, unbounded)
g_r(t)      = P(R=1|past)              # from fit_observation (raw, unbounded)
g_product(t) = g_a * g_c * g_r         # combined point-in-time probability
g_bounded(t) = bound(g_product, bounds) # bounded to prevent extreme weights
g_ac(t)     = g_bounded(t) / g_r(t)    # bounded AC component (absorbs adjustment)
sw_ac(t)    = (marg_a * marg_c) / g_ac  # stabilized AC weight
csw_ac(t)   = cumprod(sw_ac) over time  # CUMULATED (absorbing)
sw_r(t)     = marg_r / g_r(t)          # observation (NOT cumulated)
final(t)    = csw_ac(t) * sw_r(t)      # total weight
```

## TMLE algorithm (critical logic)

```
For each target time T:
1. Scale Y to [0,1] for continuous (all types use quasibinomial)
2. Compute g_cum(s) = cumprod(g_a * g_c) per subject (via .compute_cumulative_g)
3. Backward from T to min(time), at each step s:
   a. Fit Q model (quasibinomial) on risk set with non-NA Q
   b. Predict with A=regime for all risk set -> Q_bar
   c. Fluctuate: quasibinomial GLM, offset=logit(Q_bar), weights=1/(g_cum * g_r)
   d. Q*_s = expit(logit(Q_bar) + epsilon) for ALL risk set
   e. Propagate Q*_s backward as pseudo-outcome
4. psi_hat = mean(Q*_0), back-transformed if continuous
5. EIF: D_i = (Q*_0 - psi) + sum_s H_s * (Q*_{s+1} - Q*_s)
```

## Contrasts

`contrast()` is a standalone function that computes treatment effect contrasts
between two regimes. It returns a `longy_contrast` S3 object with
`print`/`summary`/`plot` methods.

```r
contrast(obj, regime = c("always", "never"))
contrast(obj, regime = "always", ref = "never")  # lmtp-style ref argument
contrast(obj, regime = c("always", "never"), scale = "ratio", estimator = "tmle")
```

### Signature

`contrast(obj, regime, ref = NULL, estimator = NULL, scale = c("difference", "ratio", "odds_ratio"), ci_level = 0.95)`

- `regime`: length-2 character vector (regime1 vs regime2), or length-1 when `ref` is provided
- `ref`: reference regime name (alternative to passing two names in `regime`)
- `estimator`: auto-detected if NULL (prefers tmle > ipw > gcomp > unadjusted)
- `scale`: `"difference"` (default), `"ratio"`, or `"odds_ratio"`

### Per-subject influence curves

Delta-method inference requires per-subject ICs stored on each `longy_result$ic`
(a `data.table` with columns: id, `.time`, `.ic`):

- **IPW**: ICs computed during `estimate_ipw(inference = "ic")` — full-population
  vector with 0 for subjects not at risk. Formula:
  `IC_i = N * w_i * (Y_i - psi) / sum(w)`
- **TMLE**: EIF values stored during `estimate_tmle(inference = "eif")` —
  one value per subject per target time
- **G-comp**: No ICs available (bootstrap only) — `contrast()` returns point
  estimates with a message, SEs are NA

### Delta method by scale

- **Difference**: `IC_contrast = IC_1 - IC_0`, `SE = sd(IC_contrast) / sqrt(n)`
- **Ratio**: `IC_ratio = (IC_1 * psi_0 - IC_0 * psi_1) / psi_0^2`
- **Odds ratio**: Delta method on log-odds scale, CIs exponentiated

### Auto-contrast via longy()

`longy(..., contrast = TRUE)` auto-computes pairwise contrasts for all regime
pairs and estimators. These are stored in `obj$contrasts` (a list keyed as
`{r1}_vs_{r2}_{estimator}`). Note: the `contrasts` field is only created by
`longy()` — it is not part of the `longy_data()` constructor.

## Cross-fitting

- Enabled via `longy(..., cross_fit = 5)` or `set_crossfit(obj, n_folds = 5)`
- Folds assigned at subject level (`.longy_fold` column)
- Dispatch at `fit_*` level: `.cf_fit_treatment()`, `.cf_fit_censoring()`,
  `.cf_fit_observation()` in `crossfit.R`
- Marginal rates computed from full risk set (population constants)
- CV-TMLE via `.cf_estimate_tmle()` with pooled fluctuation
- `.remove_tracking_columns()` preserves `.longy_fold`

## Open design decisions

- **Continuous outcome scaling**: longy currently scales continuous Y to [0,1]
  and uses quasibinomial throughout the TMLE backward pass (both Q-model fitting
  and fluctuation). This matches ltmle/stremr's approach. An alternative is
  lmtp's hybrid approach (scale to [0,1] but use gaussian for Q-models, binomial
  only for the fluctuation step). See `notes/continuous_outcome_scaling.md` for
  full comparison.

## TODO

- **Continuous outcome scaling**: Revisit whether the current approach
  (scale + quasibinomial throughout) is optimal, or whether lmtp's hybrid
  approach (scale + gaussian Q-models, binomial fluctuation only) would be
  preferable. See Open design decisions above.

## See also

- DESIGN.md for architecture decisions
- DEVLOG.md for full development history and changelog
- `ipw run.R` in parent directory for the reference implementation
