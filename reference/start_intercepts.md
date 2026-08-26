# Start at the Intercept-Only Fit

The default: each equation's intercept at the maximum likelihood
estimate of the model with every covariate removed, and every other
coefficient at zero.

## Usage

``` r
start_intercepts()
```

## Value

A `StartIntercepts` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

The model with no covariates is the model with every slope set to zero,
so its maximum likelihood estimate is exactly where the model with
covariates should begin. It costs one small fit, on an intercept-only
design.

Every other coefficient starts at zero, which is where a penalized
block's penalty is smallest.

Two adjustments make it usable where a naive intercept would not be:

- An **offset** is subtracted before the value is written. Without that,
  a model with an offset averaging 6.74 begins at \\\exp(0.897 + 6.74) =
  2080\\ against a sample mean of 2.45. Measured on a negative binomial
  over person-years: 4.9 s and 9 scoring iterations with the correction,
  more than 25 minutes without.

- The value is written to a **parametric intercept** and to nothing
  else. A term such as
  [`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
  names the intercept of each of its own parameters `(Intercept)` too,
  and those live on that parameter's chart, not the predictor's. Writing
  a predictor-scale value onto a `phi` that rides a log link gave
  \\\exp(23.9) = 2.5 \times 10^{10}\\.

This is what
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does when `start` is `NULL`. The constructor exists so the default can
be named, compared against and passed on.

## See also

[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md),
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)

## Examples

``` r
start_intercepts()
#> <start> intercept-only fit

set.seed(1)
dd <- data.frame(x = runif(60, 0, 10))
dd$y <- 500 + 20 * dd$x + rnorm(60, sd = 5)
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)

# The intercept starts at the response's own scale, and the slope at zero.
start_at(start_intercepts(), spec, statmod_design(spec), NULL)
#> $mu
#> [1] 603.0251   0.0000
#> 
#> $sigma
#> [1] 3.979791
#> 
```
