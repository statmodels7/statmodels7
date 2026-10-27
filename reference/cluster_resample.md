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

Beside a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term the groups drawn more than once become new groups: each copy of a
cluster has its own groups, with an effect of its own, as in Stata's
`bootstrap, cluster() idcluster()`. A weight is not enough in this case,
because the prior of the effects would count a group drawn twice once.
The grouping variable of each random-effect term is relabeled in the
replicate when every one of its groups lies in a single cluster (the
groups of a term that is the cluster variable itself, or nested in it).
A grouping variable that spans several clusters keeps its labels, and
its effects are shared by the copies. The design is rebuilt on the rows
drawn: the other terms are reapplied with
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html),
so a basis keeps the knots of the fit and a factor its contrasts, and
the random-effect terms are built on the new groups with the
hyperparameters of the fit.

The individual effects of different replicates are effects of different
groups, so no remapping between replicates exists. The replicates store
the coordinates that every replicate has: the coefficients of the other
terms, and with `refit = "full"` the hyperparameters, the variance
components among them. The coordinates of the individual random effects
are missing in the stored matrix, and their rows and columns in the
variance are `NA`. A term that shares a covariance through a label, a
smooth random effect
([`modelterms7::RandomSmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomSmoothTerm.html))
and a random effect inside the subformula of another term are rejected.

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

# beside a random intercept the groups drawn again are new groups
fit_re <- statmod(y ~ x + random(~ 1 | g),
                  distributions7::gaussian1_distrib(), dd)
fit_re <- add_bootstrap(fit_re, n = 10, resample = cluster_resample(~ g))
round(sqrt(diag(vcov(fit_re, vcov_type = bootstrap_vcov(
  resample = cluster_resample(~ g)))))[1:2], 3)
#> mu:(Intercept)           mu:x 
#>          0.347          0.223 
```
