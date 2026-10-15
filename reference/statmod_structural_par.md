# The Structural Terms' Estimated Parameters

One named vector per structural term, on the scale that the term's own
documentation uses for its parameters, together with the unconstrained
values that were estimated.

## Usage

``` r
statmod_structural_par(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design the fit ran on.

## Value

A named list, one entry per structural term, each a list with
`parameter` (the named vector on the parameter scale), `unconstrained`,
`equation` (the distribution parameter in which the term sits) and
`held` (the names of the parameters that an intercept of the equation
carries). Empty when there is none.

## See also

[`statmod_fit_structural()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_structural.md)
