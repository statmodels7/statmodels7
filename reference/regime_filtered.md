# The Filtered State Probabilities of a Regime Term

Runs the forward recursion of the term at its fitted parameters:
\\\alpha_t \propto f_t \circ (\alpha\_{t-1} P)\\, started from the
stationary distribution in each group, normalized at every step. The
normalizing constants are the contributions of the observations to the
log-likelihood of the term.

## Usage

``` r
regime_filtered(r)
```

## Arguments

- r:

  An entry of the `regimes` element of
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).

## Value

A list with `prob`, a matrix with one row per observation and one column
per regime, and `loglik`, the contributions of the observations.
