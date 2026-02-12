# longy Design Document

## Estimator Roadmap

### v0.1 (current)
- IPW with stabilized weights
- IC-based inference (+ sandwich, bootstrap)
- Weight diagnostics and positivity diagnostics
- Three nuisance models: g_A, g_C, g_R

### v0.2 (planned)
- G-computation (outcome regression)
- Sequential regression / iterated conditional expectations

### v0.3 (planned)
- TMLE (targeted minimum loss-based estimation)
- Doubly-robust estimator combining IPW + G-comp
- Cross-fitting infrastructure activated

### v0.4+ (planned)
- Stochastic interventions
- Continuous treatments
- Competing risks
- Marginal structural models

## Data Model

Long format only: one row per (subject, time). All variables are columns.

### Node types
- **id**: subject identifier
- **time**: integer time index
- **outcome** (Y): the event/measurement of interest
- **treatment** (A): binary intervention variable
- **censoring** (C): absorbing dropout indicators (can be multiple). When multiple
  censoring sources exist, each is modeled independently and weights are multiplied.
  The package does not model within-interval ordering of censoring events — if
  multiple sources can trigger in the same interval, the data should record only
  the first censoring event.
- **observation** (R): intermittent outcome measurement indicator
- **baseline** (W): time-invariant covariates
- **timevarying** (L): time-varying covariates

### Censoring vs Observation
This is a critical distinction:
- **Censoring (C)**: absorbing. Once C=1, the subject is gone forever. Weights are cumulated.
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
- Uncensored at time t (C(t) = 0)

## Key Design Decisions

1. **S3 over R6**: Keeps things simple and testable. Users don't need to understand
   reference semantics.

2. **data.table internally**: Performance for large longitudinal datasets. Convert
   at boundary (accept data.frame, return data.frame in results).

3. **SuperLearner in Suggests**: Package works with just glm. SuperLearner is optional
   for users who want ML-based nuisance estimation.

4. **Pipe-friendly**: Each step returns the modified longy_data object, enabling
   `|>` chaining.

5. **Cross-fitting deferred**: Architecture supports it (fold_id in longy_data),
   but actual sample-splitting is v0.3.

6. **No ffSL dependency**: The reference code uses a custom parallel SL wrapper.
   We use standard SuperLearner with glm fallback. Users can extend later.

## Open Questions

- Should we support multiple outcome types (survival, continuous, binary) in v0.1?
  Decision: Yes, outcome_type is stored but v0.1 mainly targets binary.
- Marginal structural model (MSM) smoothing across time points?
  Decision: Deferred to v0.4.
