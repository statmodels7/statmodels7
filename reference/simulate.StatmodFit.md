# Draws from a Fitted Model

Responses drawn from the fitted distribution, one column per replicate.

## Usage

``` r
# S3 method for class 'StatmodFit'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- nsim:

  How many replicates.

- seed:

  Passed to [`base::set.seed()`](https://rdrr.io/r/base/Random.html) if
  given, the caller's stream being restored afterwards.

- ...:

  Unused.

## Value

A data frame with `nobs(object)` rows and `nsim` columns, each column
one draw of the whole response vector.

## Details

The draws are taken at the fitted parameters, so what they carry is the
variation of the response alone. The uncertainty of the estimates is not
in them; a parametric bootstrap adds it by refitting each column.

For a model carrying a structural term the fitted parameters are the
ones the filter reached along the observed series, so the draws are
conditional on that series and are not a fresh path of the process.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
to simulate from a model that was never fitted,
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the fitted parameters the draws come from

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
dim(simulate(fit, nsim = 3))
#> [1] 40  3
```
