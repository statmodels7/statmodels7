# Start at the Intercept-Only Fit

The default: each equation's intercept at the maximum likelihood
estimate of the model with every covariate removed, and every other
coefficient at zero.

## Usage

``` r
start_intercepts()
```

## Value

A
[`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

The model with no covariates is the same model with every slope set to
zero, so it is exactly where the model with them should begin, and it
costs one small fit. The penalized blocks start at zero, where their
penalty is smallest. This is what
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
does when `start` is `NULL`; the constructor exists so that the default
can be named, compared against and passed on.

## See also

[`start_origin`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random`](https://statmodels7.github.io/statmodels7/reference/start_random.md)

## Examples

``` r
start_intercepts()
#> <start> intercept-only fit
```
