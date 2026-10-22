# Run a Joint Filter at the Current Parameters

Returns the entry of
[`statmod_filter_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
for a joint filter: the filtered predictors of every filtered equation,
their Jacobian in the filter's parameters, and the callbacks that
produced them.

## Usage

``` r
joint_filter_entry(spec, u, tm, psi, eta_static)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- u:

  The structural unit.

- tm:

  The joint filter.

- psi:

  Its parameters on the parameter scale, as
  [`structural_psi()`](https://statmodels7.github.io/statmodels7/reference/structural_psi.md)
  returns them.

- eta_static:

  The static predictors, a named list.

## Value

A list with `param` (the first filtered parameter), `params`, `joint`
(`TRUE`), `term`, `tm`, `eta` (a matrix with one column per filtered
parameter), `jacobian` (a named list of matrices), `curv`, `psi`, `cb`
and `eta_static` (a matrix like `eta`).
