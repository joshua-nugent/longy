# longy — Claude Code Context

## What is longy?

longy ("lawn-gee") is a private R package for longitudinal causal inference — an
alternative to `ltmle` and `stremr`. It provides IPW, G-computation, and TMLE
estimators with IC/EIF/bootstrap/sandwich inference, plus cross-fitting (CV-TMLE).

## Architecture

- **S3 classes**: `longy_data` (single accumulating pipeline object), `longy_result` (lightweight result)
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
                    observation, sampling_weights,
                    baseline, timevarying, outcome_type, competing),
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
| R/estimate_ipw.R | `estimate_ipw()` + print/summary/plot for all estimators |
| R/estimate_gcomp.R | `estimate_gcomp()` — G-comp estimator + bootstrap |
| R/estimate_tmle.R | `estimate_tmle()` — TMLE estimator + EIF inference |
| R/inference.R | IC-based + bootstrap inference (IPW, G-comp, TMLE) |
| R/diagnostics.R | Weight & positivity diagnostics |
| R/crossfit.R | Cross-fitting: `.cf_fit_treatment/censoring/observation()`, `.cf_estimate_tmle()` |
| R/longy.R | High-level `longy()` wrapper + `results()` accessor + `plot.longy_data` |
| R/ffSL.R | Future-factorial SuperLearner (parallel CV via `future.apply`) |
| R/utils.R | Internal helpers: `.as_longy_data()`, `.safe_sl()`, `.bound()`, `.predict_from_fit()`, `.vmsg()` |
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

```
sw_a(t) = marginal_a / p_a(t)          # treatment
sw_c(t) = marginal_c / p_c(t)          # censoring (per cause)
sw_ac(t) = sw_a(t) * sw_c(t)           # combined point-in-time
csw_ac(t) = cumprod(sw_ac) over time   # CUMULATED (absorbing)
sw_r(t) = marginal_r / p_r(t)          # observation (NOT cumulated)
final(t) = csw_ac(t) * sw_r(t)         # total weight
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

## Cross-fitting

- Enabled via `longy(..., cross_fit = 5)` or `set_crossfit(obj, n_folds = 5)`
- Folds assigned at subject level (`.longy_fold` column)
- Dispatch at `fit_*` level: `.cf_fit_treatment()`, `.cf_fit_censoring()`,
  `.cf_fit_observation()` in `crossfit.R`
- Marginal rates computed from full risk set (population constants)
- CV-TMLE via `.cf_estimate_tmle()` with pooled fluctuation
- `.remove_tracking_columns()` preserves `.longy_fold`

## Open design decisions

- **Continuous outcome scaling**: longy currently uses gaussian (no scaling) for
  continuous outcomes. lmtp scales to [0,1] but still uses gaussian for Q-models
  (binomial only for TMLE fluctuation). ltmle/stremr scale to [0,1] with
  quasibinomial throughout. Revisit whether to adopt lmtp's hybrid approach
  (scale + gaussian). See `notes/continuous_outcome_scaling.md` for full comparison.

## TODO

- **Quasibinomial learner compatibility**: Backward ICE pseudo-outcomes are
  continuous [0,1], requiring quasibinomial family. Many learners (glmnet,
  xgboost binary:logistic, others) crash or misbehave with non-integer Y.
  Current workarounds are per-learner (swap xgboost to reg:squarederror, drop
  glmnet when Y is non-binary). Need a more general solution — possibly a
  wrapper that intercepts any binomial-family learner and substitutes a
  gaussian-family version with truncated predictions.
- **Continuous outcome scaling**: Decide whether to scale Y to [0,1] for
  continuous outcomes (see Open design decisions above).

## See also

- DESIGN.md for architecture decisions
- DEVLOG.md for full development history and changelog
- `ipw run.R` in parent directory for the reference implementation
