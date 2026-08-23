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
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- nsim:

  How many replicates.

- seed:

  Passed to [`set.seed`](https://rdrr.io/r/base/Random.html) if given,
  the caller's stream being restored afterwards.

- ...:

  Unused.

## Value

A data frame of `nsim` columns.

## Details

The draws are taken at the FITTED parameters, so what they carry is the
variation of the response and not the uncertainty of the estimates; a
parametric bootstrap adds the second by refitting each column.

For a model carrying a structural term the fitted parameters are the
ones the filter reached ALONG the observed series, so the draws are
conditional on that series rather than a fresh path of the process.

## See also

[`rstatmod`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md),
[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
dim(simulate(fit, nsim = 3))
#> [1] 40  3
```
