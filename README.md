# longy

Longitudinal causal inference in R. Estimates causal effects with time-varying treatments, confounders, informative censoring, and intermittent outcome measurement.

Supports IPW (inverse probability weighting) and G-computation estimators. TMLE coming soon.

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
  estimator = "ipw",
  n_boot = 0,
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

## G-computation

```r
obj <- longy_data(
  sim_longy,
  id = "id", time = "time", outcome = "Y",
  treatment = "A", censoring = "C", observation = "R",
  baseline = c("W1", "W2"), timevarying = c("L1", "L2")
)

obj <- define_regime(obj, "always", static = 1L)
obj <- fit_outcome(obj, regime = "always")
result <- estimate_gcomp(obj, regime = "always", n_boot = 0)

result
```

## Features

- **Long format** data (one row per person-time)
- **IPW** with stabilized weights, influence-curve / bootstrap / sandwich inference
- **G-computation** via iterated conditional expectations (backward sequential regression)
- **Binary, continuous, and survival** outcomes
- **Informative censoring** and **intermittent missingness** handling
- **Static and dynamic** treatment regimes
- Optional **SuperLearner** ensemble learning
- Parallel bootstrap via **future.apply**
