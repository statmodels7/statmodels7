# Values of a Profile

Selects the values at which
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
evaluates the profile of a coefficient: `n` equally spaced values from
`width` standard errors below the estimate to `width` standard errors
above it, on the scale of the coefficient in the linear predictor.

## Usage

``` r
profile_range(width = 3, n = 21)

ProfileRange(width = integer(0), n = integer(0))
```

## Arguments

- width:

  A single positive number: the half-width of the grid in standard
  errors.

- n:

  A single integer, at least 2: the number of values.

## Value

An object of class `ProfileRange`, with the properties `width` and `n`.

## Details

The standard error is the one of `vcov(fit)`, so it follows the variance
convention that [`vcov()`](https://rdrr.io/r/stats/vcov.html) uses by
default. The grid is symmetric about the estimate, and it contains the
estimate itself when `n` is odd.

The scale is the unconstrained one. A coefficient of a linear predictor
lives there, and so do the parameters of a nonlinear term, the positions
of the break-points and the parameters of a covariance structure; the
scale of the coefficients is described at
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md).

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
which takes it in its argument `values`.

## Examples

``` r
profile_range()
#> <statmodels7::ProfileRange>
#>  @ width: num 3
#>  @ n    : int 21
profile_range(width = 2, n = 9)
#> <statmodels7::ProfileRange>
#>  @ width: num 2
#>  @ n    : int 9
```
