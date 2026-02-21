# longy Design Document

## Estimator Roadmap

### v0.1 (completed 2026-02-12)
- IPW with stabilized weights
- IC-based inference (+ sandwich, bootstrap)
- Weight diagnostics and positivity diagnostics
- Three nuisance models: g_A, g_C, g_R

### v0.2 (completed)
- G-computation (outcome regression)
- Sequential regression / iterated conditional expectations (backward ICE)

### v0.3 (completed)
- TMLE (targeted minimum loss-based estimation)
- Doubly-robust estimator via backward ICE + quasibinomial fluctuation + EIF inference
- Competing risks support (`competing` parameter)

### v0.4 (completed 2026-02-18)
- Cross-fitting infrastructure (CV-TMLE with pooled fluctuation)
- g_R in TMLE clever covariate: H(s) = 1 / (g_cum(s) * g_r(s))
- ffSL (future-factorial parallel SuperLearner)

### v0.5 (completed 2026-02-21)
- Single character column censoring interface (replaces multiple binary columns)
- Multi-cause censoring decomposition (`"death"`, `"ltfu"`, etc.)

### Future
- Stochastic interventions
- Continuous treatments
- Marginal structural models

## Data Model

Long format only: one row per (subject, time). All variables are columns.

### Node types
- **id**: subject identifier
- **time**: integer time index
- **outcome** (Y): the event/measurement of interest
- **treatment** (A): binary intervention variable (0/1)
- **censoring** (C): a single character/factor column with values like `"uncensored"`,
  `"censored"`, `"death"`, `"ltfu"`. Internally decomposed into binary `.cens_<cause>`
  indicator columns by `longy_data()`. Each cause is modeled independently and
  their weights are multiplied. The package does not model within-interval
  ordering of censoring events.
- **observation** (R): intermittent outcome measurement indicator (binary 0/1)
- **baseline** (W): time-invariant covariates
- **timevarying** (L): time-varying covariates
- **competing** (D): binary absorbing competing event indicator (0/1), used with
  `outcome_type = "survival"` only
- **sampling_weights**: external survey/sampling weights (numeric, constant within subject)

### Censoring vs Observation
This is a critical distinction:
- **Censoring (C)**: absorbing. Once a non-`"uncensored"` value appears, the subject is
  gone forever. Weights are cumulated over time.
- **Observation (R)**: intermittent. R=0 means outcome not measured at this time, but
  subject can return at t+1. Weights are point-in-time (NOT cumulated).

## Risk Set Logic

Each nuisance model has a specific risk set:

### g_A (treatment)
At-risk at time t if:
- Regime-consistent through t-1 (followed the prescribed treatment through previous time)
- Uncensored through t-1

### g_C (censoring)
At-risk at time t if:
- Regime-consistent through t-1
- Uncensored through t-1
- Treatment at time t consistent with regime (already received A(t))

### g_R (observation)
At-risk at time t if (most restrictive):
- Regime-consistent through t-1
- Uncensored through t-1
- Treatment at time t consistent with regime
- Uncensored at time t (C(t) = "uncensored")

## Key Design Decisions

1. **S3 over R6**: Keeps things simple and testable. Users don't need to understand
   reference semantics.

2. **data.table internally**: Performance for large longitudinal datasets. Convert
   at boundary (accept data.frame, return data.frame in results).

3. **SuperLearner in Suggests**: Package works with just glm. SuperLearner is optional
   for users who want ML-based nuisance estimation.

4. **Pipe-friendly**: Each step returns the modified longy_data object, enabling
   `|>` chaining.

5. **Cross-fitting at subject level**: Folds are assigned to subjects, not rows.
   All rows for a subject belong to the same fold. Implemented via `.longy_fold`
   column and dispatch at the `fit_*` level.

6. **ffSL for parallelism**: The `sl_fn = "ffSL"` option parallelizes SuperLearner
   CV folds x algorithms via `future.apply`. Default in `longy()` wrapper;
   individual `fit_*` functions default to `"SuperLearner"` (sequential).

7. **Character censoring column**: Users provide a single character column
   (e.g. `"uncensored"`, `"death"`, `"ltfu"`). `longy_data()` decomposes it into
   binary `.cens_<cause>` columns internally. This keeps the user interface simple
   while all downstream code operates on binary indicators as before.

## Outcome Types

Three outcome types are supported:
- **binary**: Y in {0, 1}
- **continuous**: Y numeric, scaled to [0,1] internally for TMLE (quasibinomial)
- **survival**: Y absorbing (once Y=1, stays 1), optionally with competing risks
