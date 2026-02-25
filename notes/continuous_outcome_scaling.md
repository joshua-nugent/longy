# Continuous Outcome Handling: Package Comparison

## Date: 2026-02-24

## Context

When the outcome Y is continuous, packages must decide:
1. Whether to scale Y to [0,1] before modeling
2. What family/link to use for Q-model regressions (outcome + backward ICE pseudo-outcomes)
3. How to bound predictions to prevent extrapolation

## Package approaches

### ltmle
- **Scales Y to [0,1]** using observed `Yrange = range(Y)`
- **quasibinomial (logit link)** for ALL Q-models (outcome + pseudo-outcomes)
- Restricts learner compatibility (glmnet, xgboost crash on continuous [0,1] with binomial)

### stremr
- **Requires Y pre-scaled to [0,1]** by the user
- **quasibinomial (logit link)** throughout
- `maxpY` parameter for rescaling (default 1.0)
- Enforces strict [0,1] bounds at every backward step
- Same learner restrictions as ltmle

### tmle (classic, single time-point)
- **Scales Y to [0,1]** using `range(Y)` when `fluctuation="logistic"` (default)
- Initial Q-model fit with **gaussian**, then transforms to logit scale
- Fluctuation step uses **binomial (logit link)**
- Back-transforms to original scale at the end

### lmtp (modern)
- **Scales Y to [0,1]** using observed range
- **gaussian family** for Q-model regressions (full learner compatibility)
- At intermediate backward steps, outcome type forced to "continuous" → gaussian
- Only the **TMLE fluctuation step** uses binomial (logit link) with `qlogis(offset)`
- Predictions bounded to [1e-5, 0.99999]

### longy (current)
- **No scaling** — models Y on the raw scale
- **gaussian (identity link)** for all Q-models
- No bounding of predictions
- Full learner compatibility

## Summary table

| Package | Y scaling | Q-model family | Fluctuation family | Learner compat |
|---------|-----------|----------------|-------------------|----------------|
| ltmle   | [0,1]     | quasibinomial  | quasibinomial     | limited        |
| stremr  | [0,1]     | quasibinomial  | quasibinomial     | limited        |
| tmle    | [0,1]     | gaussian→logit | binomial          | mixed          |
| lmtp    | [0,1]     | gaussian       | binomial          | full           |
| longy   | none      | gaussian       | N/A yet           | full           |

## Empirical finding (debug_gcomp_continuous.R)

In a simcausal comparison with `Y ~ N(mean, sd=1)`:
- longy (gaussian, no scaling) was **closer to MC truth** than ltmle (quasibinomial, scaled)
- The logit link introduces some distortion for truly continuous outcomes

## Recommendation (deferred)

Consider adopting lmtp's hybrid approach:
1. Scale Y to [0,1] using observed range (bounds pseudo-outcomes, prevents extrapolation)
2. Keep gaussian for Q-models (preserves full learner compatibility)
3. Use binomial for TMLE fluctuation step (theoretically motivated)
4. Back-transform estimates to original scale

This gives bounded predictions without sacrificing learner compatibility. Decision deferred —
current gaussian-no-scaling approach works well in testing.
