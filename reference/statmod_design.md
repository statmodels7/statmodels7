# The Design of a Specification

Returns, per parameter, the terms' blocks side by side and the names of
the coefficients they carry.

## Usage

``` r
statmod_design(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A named list with one entry per parameter, each a list with `X`,
`coef_names`, `npar` and `blocks` (the column range each term occupies).
Where no term moves with the coefficients and none is structural, the
attribute `eta_memo` is an environment holding the last predictors
[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
computed on this design, so a second call at the same coefficients
returns them without recomputing.

## See also

[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)

## Examples

``` r
dd <- data.frame(y = rnorm(20), x = runif(20))
d <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
vapply(statmod_design(d), function(z) z$npar, integer(1))
#>    mu sigma 
#>     2     1 
```
