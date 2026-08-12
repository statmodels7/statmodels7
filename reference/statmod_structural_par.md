# The Structural Terms' Estimated Parameters

One named vector per structural term, on the scale the term's own
documentation describes them on, together with the unconstrained values
they were estimated as.

## Usage

``` r
statmod_structural_par(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design the fit ran on.

## Value

A named list, one entry per structural term, each a list with
`parameter` and `unconstrained`. Empty when there is none.

## See also

[`statmod_fit_structural`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md)
