# Read the Stored Bootstrap of a Fit

Returns the table of the bootstrap entries that
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
stored in a fit, or one entry.

## Usage

``` r
statmod_bootstrap(fit, resample = NULL, refit = NULL)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  A
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
  object, or `NULL`.

- refit:

  `"coefficients"`, `"full"` or `NULL`.

## Value

A data frame, or a
[StatmodBootstrap](https://statmodels7.github.io/statmodels7/reference/StatmodBootstrap.md)
object.

## Details

Without `resample` and `refit` the result is a data frame with one row
for each entry: the resampling scheme, the refit, the number of
replicates, and the number of replaced replicates in all and by cause.
With them the result is the entry that has those settings, a
[StatmodBootstrap](https://statmodels7.github.io/statmodels7/reference/StatmodBootstrap.md)
object, and an argument that is `NULL` matches any value. If no entry
matches, or if more than one does, an error is signaled; the second one
prints the table. [`coef()`](https://rdrr.io/r/stats/coef.html) of an
entry is the matrix of replicates.

## See also

[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md),
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
statmod_bootstrap(fit)
#> [1] resample     refit        replicates   replaced     error       
#> [6] mode_error   unidentified
#> <0 rows> (or 0-length row.names)
fit <- add_bootstrap(fit, n = 20, resample = case_resample())
statmod_bootstrap(fit)
#>   resample        refit replicates replaced error mode_error unidentified
#> 1     case coefficients         20        0     0          0            0
head(coef(statmod_bootstrap(fit, resample = case_resample())), 3)
#>      mu:(Intercept)     mu:x sigma:(Intercept)
#> [1,]       1.135981 1.753653        -0.9760249
#> [2,]       1.242773 1.771733        -1.0491633
#> [3,]       1.207261 1.753075        -0.9625627
```
