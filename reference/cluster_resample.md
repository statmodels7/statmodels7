# Cluster Resampling

Selects the cluster bootstrap: every replicate draws the groups of
observations with replacement and refits the model on them.

## Usage

``` r
cluster_resample(cluster, data = NULL)
```

## Arguments

- cluster:

  A one-sided formula with one variable, as in `~ id`.

- data:

  A data frame with that variable and one row for each observation of
  the fit, or `NULL`. With `NULL` the data frame is the `data` of the
  call that built the fit.

## Value

A `ClusterResample` object, inheriting from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md),
with the properties `cluster` and `data`.

## Details

The groups are the values of the variable that the formula names, a
column of `data`. A fit does not keep its data: with `data = NULL` the
data frame is the `data` argument of the call of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
evaluated in the environment of the model formula, and an error says to
pass `data` where that fails or where its number of rows is not the
fit's. The variable must have no missing values. A group drawn \\k\\
times gives all its observations the count \\k\\. Without random effects
the likelihood is a sum over observations, so this is the fit to the
duplicated groups with the design held.

In this version a fit with a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term is rejected. The coordinates of the individual random effects have
no counterpart between replicates, because a group can be absent or
appear twice.

## See also

[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
for the comparison of the three schemes,
[`add_bootstrap()`](https://statmodels7.github.io/statmodels7/reference/add_bootstrap.md)
and
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
for the functions that take it.

## Examples

``` r
set.seed(1)
dd <- data.frame(g = rep(1:12, each = 5), x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(12)[dd$g] + rnorm(60, sd = 0.4)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit <- add_bootstrap(fit, n = 30, resample = cluster_resample(~ g))
statmod_bootstrap(fit)
#>       resample        refit replicates replaced error mode_error unidentified
#> 1 cluster (~g) coefficients         30        0     0          0            0
```
