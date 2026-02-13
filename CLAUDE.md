# longy — Claude Code Context

## What is longy?

longy ("lawn-gee") is a private R package for longitudinal causal inference — an
alternative to `ltmle` and `stremr`. v0.1 provides IPW estimation with IC-based
inference; future versions will add G-computation, TMLE, and cross-fitting.

## Architecture

- **S3 classes** for user-facing objects (`longy_data`, `longy_result`)
- **Plain functions** internally (no R6)
- **data.table** internally for performance; accept data.frame at boundary
- **SuperLearner** for ML with glm fallback (SuperLearner is in Suggests, not Imports)
- **Long format only**: one row per person-time

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

## Key data structure: longy_data object

```r
list(
  data       = <data.table keyed on (id, time)>,
  nodes      = list(id, time, outcome, treatment, censoring, observation,
                    baseline, timevarying, outcome_type, competing_risks),
  regimes    = list(),
  fits       = list(treatment = NULL, censoring = list(), observation = NULL),
  weights    = NULL,
  crossfit   = list(enabled = FALSE, n_folds = NULL, fold_id = NULL),
  meta       = list(n_subjects, n_obs, n_times, time_values, max_time, min_time)
)
```

## File map

| File | Purpose |
|------|---------|
| R/data_input.R | `longy_data()` constructor + validation |
| R/regimes.R | `define_regime()` + `.evaluate_regime()` |
| R/fit_treatment.R | g_A models |
| R/fit_censoring.R | g_C models |
| R/fit_observation.R | g_R models (intermittent) |
| R/weights.R | `compute_weights()` |
| R/fit_outcome.R | Outcome models via sequential regression (G-comp) |
| R/estimate_ipw.R | `estimate_ipw()` |
| R/estimate_gcomp.R | `estimate_gcomp()` — G-comp estimator + bootstrap |
| R/inference.R | IC-based + bootstrap inference |
| R/diagnostics.R | Weight & positivity diagnostics |
| R/crossfit.R | Cross-fitting infrastructure (stub in v0.1) |
| R/longy.R | High-level `longy()` wrapper |
| R/utils.R | Internal helpers |

## Coding conventions

- snake_case for everything
- Internal helpers prefixed with `.` (e.g., `.bound()`, `.safe_sl()`)
- data.table internally, data.frame at boundaries
- All exported functions get roxygen2 docs
- Tests use testthat edition 3
- SuperLearner access via `requireNamespace("SuperLearner", quietly = TRUE)`

## Weight computation (critical logic)

```
sw_a(t) = marginal_a / p_a(t)          # treatment
sw_c(t) = marginal_c / p_c(t)          # censoring (per source)
sw_ac(t) = sw_a(t) * sw_c(t)           # combined point-in-time
csw_ac(t) = cumprod(sw_ac) over time   # CUMULATED (absorbing)
sw_r(t) = marginal_r / p_r(t)          # observation (NOT cumulated)
final(t) = csw_ac(t) * sw_r(t)         # total weight
```

## See also

- DESIGN.md for full architecture decisions and estimator roadmap
- `ipw run.R` in parent directory for the reference implementation
