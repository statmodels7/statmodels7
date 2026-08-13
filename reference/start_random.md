# Start From a Random Draw

Every coefficient drawn on the unconstrained scale, by default added to
the intercept-only fit rather than replacing it.

## Usage

``` r
start_random(fn = stats::rnorm, ..., center = TRUE)
```

## Arguments

- fn:

  A generator called as `fn(n, ...)`, returning `n` values. Defaults to
  [`rnorm`](https://rdrr.io/r/stats/Normal.html).

- ...:

  Further arguments to `fn`, such as `sd` or `min` and `max`.

- center:

  Whether to add the draw to the intercept-only start rather than use it
  alone. Defaults to `TRUE`.

## Value

A
[`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

The draw is added to
[`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)'s
answer unless `center = FALSE`, and that is not a detail: a coefficient
drawn from a standard normal is a sensible perturbation and a hopeless
absolute value, since the intercept of a location equation is on the
scale of the response. Centring keeps the scale and randomizes the
direction, which is what a caller wanting several starts is after.

The stream is the caller's, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the result
and a fit begun this way is reproducible only alongside its seed.

## See also

[`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.html)

## Examples

``` r
start_random()
#> <start> random around the intercept-only fit
start_random(stats::runif, min = -2, max = 2)
#> <start> random around the intercept-only fit
```
