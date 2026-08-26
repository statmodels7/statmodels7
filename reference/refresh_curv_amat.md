# The Weight the Second Derivative of a Block Is Paired With

[`refresh_amat()`](https://statmodels7.github.io/statmodels7/reference/refresh_amat.md)
with the curvature, one matrix per unit.

## Usage

``` r
refresh_curv_amat(spec, design, M, params, npar, offs, units, Hl)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- M:

  The matrix the trace is taken against.

- params, npar, offs:

  The block bookkeeping.

- units:

  The refreshable terms, from
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- Hl:

  The link-scale curvature, from
  [`refresh_hessian()`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md).

## Value

A list of matrices, one per unit.

## Details

It is the same weight
[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
builds for the gradient, and it depends on neither direction, so it is
built once per Hessian rather than once per pair of hyperparameters.

## See also

[`trace_refresh4()`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md),
[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
