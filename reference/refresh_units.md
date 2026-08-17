# The Terms Whose Block Moves With Their Coefficients

One entry per refreshable term, carrying everything its derivative
needs: the distribution parameter it sits in, its columns in the stacked
coefficient vector, its own coefficients and the term rebuilt at them.

## Usage

``` r
refresh_units(spec, design, coef, params, npar, offs)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at `coef`.

- coef:

  The coefficients.

- params, npar, offs:

  The block bookkeeping.

## Value

A list, empty where no block moves.

## Details

Resolved once and passed to every consumer, so the gradient's
contraction and the Hessian's three corrections cannot disagree about
which terms move or where their columns are.

## See also

[`u_refresh`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
[`contract3_refresh`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md)
