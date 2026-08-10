# Fit the Smooth Block

Runs `inner_method` on the jointly fitted coefficients, the others held
fixed.

## Usage

``` r
fit_smooth(obj, beta, idx, spec, design, hyper, method, vb)
```

## Arguments

- obj:

  The objective.

- beta:

  The current stacked coefficients.

- idx:

  The smooth block's indices.

- spec:

  The specification.

- design:

  The design.

- hyper:

  The hyperparameters.

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizer.

- vb:

  The resolved verbosity.

## Value

A list with `par`, `value`, `converged`, `iterations` and `history`.
