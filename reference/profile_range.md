# Values of a Profile

Selects the values at which
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
evaluates the profile of a coefficient: `n` equally spaced values from
`width` standard errors below the estimate to `width` standard errors
above it, on the scale of the coefficient in the linear predictor.

## Usage

``` r
profile_range(width = 3, n = 21, fallback = 3)

ProfileRange(width = integer(0), n = integer(0), fallback = 3)
```

## Arguments

- width:

  A single positive number: the half-width of the grid in standard
  errors.

- n:

  A single integer, at least 2: the number of values.

- fallback:

  A single positive number: the half-width of the grid of a
  hyperparameter on the scale of its link, used where its standard error
  is not available.

## Value

An object of class `ProfileRange`, with the properties `width`, `n` and
`fallback`.

## Details

The standard error is the one of `vcov(fit)`, so it follows the variance
convention that [`vcov()`](https://rdrr.io/r/stats/vcov.html) uses by
default. The grid is symmetric about the estimate, and it contains the
estimate itself when `n` is odd.

For a hyperparameter the grid is placed on the scale of its link (for
example \\\log\lambda\\), at the estimate plus and minus `width`
standard errors, with the standard error from the exact Hessian of the
outer criterion, the one that
[`summary()`](https://rdrr.io/r/base/summary.html) reads. Where that
standard error is not available, or the estimate lies outside the
interval below, because the hyperparameter is at an edge of its range or
the criterion is flat in its direction, the grid is the estimate, cut to
that interval, plus and minus `fallback`. The values are cut to the
interval from -8 to 8 of the scale of the link, which
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads as the interior.

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
#>  @ width   : num 3
#>  @ n       : int 21
#>  @ fallback: num 3
profile_range(width = 2, n = 9)
#> <statmodels7::ProfileRange>
#>  @ width   : num 2
#>  @ n       : int 9
#>  @ fallback: num 3
```
