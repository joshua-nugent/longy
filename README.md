# longy

Longitudinal causal inference in R. Estimates causal effects with time-varying treatments, confounders, informative censoring, and intermittent outcome measurement.

Implements three estimators: **IPW** (inverse probability weighting), **G-computation**, and **TMLE** (targeted minimum loss-based estimation).

## Installation

```r
# install.packages("remotes")
remotes::install_github("joshua-nugent/longy")
```

## Quick start

```r
library(longy)

# Use the built-in simulated dataset
data(sim_longy)

# One-line analysis with the longy() wrapper
results <- longy(
  data = sim_longy,
  id = "id", time = "time", outcome = "Y",
  treatment = "A", censoring = "C", observation = "R",
  baseline = c("W1", "W2"), timevarying = c("L1", "L2"),
  regimes = list(always = 1L, never = 0L),
  estimator = "tmle",
  verbose = FALSE
)

results
```

## Modular pipeline

For more control, use the step-by-step API:

```r
obj <- longy_data(
  sim_longy,
  id = "id", time = "time", outcome = "Y",
  treatment = "A", censoring = "C", observation = "R",
  baseline = c("W1", "W2"), timevarying = c("L1", "L2")
)

obj <- define_regime(obj, "always", static = 1L)
obj <- fit_treatment(obj, regime = "always")
obj <- fit_censoring(obj, regime = "always")
obj <- fit_observation(obj, regime = "always")
obj <- compute_weights(obj, regime = "always")
result <- estimate_ipw(obj, regime = "always")

result
```

## Estimators

### IPW

Inverse probability weighting with stabilized weights. Inference via influence curves, bootstrap, or sandwich (survey) standard errors.

```r
results <- longy(sim_longy, ..., estimator = "ipw")
```

### G-computation

Outcome regression via iterated conditional expectations (backward sequential regression). Inference via bootstrap.

```r
results <- longy(sim_longy, ..., estimator = "gcomp")
```

### TMLE

Doubly-robust targeted minimum loss-based estimation. Combines outcome regression with a fluctuation step using propensity scores. Consistent if either the outcome or treatment/censoring models are correctly specified. Inference via the efficient influence function (EIF).

```r
results <- longy(sim_longy, ..., estimator = "tmle")
```

### Run multiple estimators

```r
results <- longy(sim_longy, ..., estimator = "both")  # IPW + TMLE
```

## Features

- **Long format** data (one row per person-time)
- **Three estimators**: IPW, G-computation, and TMLE
- **Binary, continuous, and survival** outcomes
- **Informative censoring** and **intermittent missingness** handling
- **Static and dynamic** treatment regimes
- **Doubly-robust** inference via TMLE with EIF-based standard errors
- Optional **SuperLearner** ensemble learning
- Parallel bootstrap via **future.apply**
- Weight and positivity **diagnostics**

## Status

v0.3 — 270 tests passing, 0 errors / 0 warnings / 0 notes on R CMD check. Cross-fitting (sample-split nuisance estimation) is next.
