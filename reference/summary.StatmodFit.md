# Summarize a Fitted Model

One coefficient table per distribution parameter – estimate, standard
error, Wald statistic, p-value and interval – with the degrees of
freedom, the information criteria and the qualifications the numbers
carry.

## Usage

``` r
# S3 method for class 'StatmodFit'
summary(
  object,
  level = 0.95,
  type = c("bayesian", "frequentist"),
  correct = FALSE,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- level:

  The confidence level.

- type:

  Which variance matrix: passed to
  [`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

- correct:

  Whether the degrees of freedom carry what the estimation of the
  hyperparameters cost. The ordinary count reads them as known, and they
  were chosen from the same data, so a criterion built on it is too
  generous. See
  [`statmod_edf_correction`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md).
  Defaults to `FALSE` because it changes a number a reader may be
  comparing with an earlier fit; it is zero where no hyperparameter was
  estimated.

- ...:

  Passed to
  [`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A
[`StatmodSummary`](https://statmodels7.github.io/statmodels7/reference/StatmodSummary-class.md).

## Details

**Each distribution parameter is read as blocks, not as one list of
coefficients**, because most of a fitted model's coefficients are not
quantities anybody reads. The blocks are

- the parametric terms:

  every unpenalized term together, one row per coefficient, which is the
  ordinary table.

- one block per smooth:

  the linear component's coefficient where the construction carries one,
  the smoothing parameter, and the effective degrees of freedom. The
  coefficients of the wiggly part are not shown: they are coordinates in
  an orthonormal basis and say nothing one at a time, while what they
  say together is the edf.

- one block per random effect:

  the parameters of the effects' distribution – what is usually called
  the variance component – and the edf. Not the effects themselves, of
  which there is one per level.

- one block per selection:

  a lasso, a SCAD or an MCP: its hyperparameters, how many coefficients
  survived, and those coefficients. The ones set exactly to zero are
  counted, not listed.

- one block per other penalized term:

  its coefficients, which stay interpretable under a ridge, together
  with its hyperparameters.

**A hyperparameter is the first row of its block**, since it governs
every coefficient under it, and the cell where its standard error would
be says what put the value there. One estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
maximizes a twice differentiable criterion, so it carries a standard
error and an interval, both read on the free scale its link defines and
mapped back
([`statmod_hyper_vcov`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_vcov.md)).
One chosen by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) over
a kinked penalty is the argument of a minimum over a grid, so the row
names the criterion and leaves the remaining columns empty: there is no
curvature at such a point to read a standard error from. One the caller
set is marked fixed.

**What a Wald p-value means here depends on the row**, and the summary
says which is which rather than printing one column and leaving it at
that. For an unpenalized coefficient it is the usual thing. For a
coefficient in a penalized block it is conditional on the smoothing
parameter, which was not estimated jointly with it, and it does not
account for the shrinkage of the estimate towards zero. For a
coefficient a kinked penalty selected, the row exists only because that
coefficient survived the selection, and a naive interval there
under-covers.

**The degrees of freedom** are the effective ones, summed over the
terms, so that a penalized term counts what it spends rather than how
many columns it has. The information criteria are built on that count.

## See also

[`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
[`confint.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(120))
dd$y <- 1 + 2 * dd$x + rnorm(120, sd = 0.4)
summary(statmod(y ~ x | sigma ~ x,
                distributions7::gaussian1_distrib(), dd))
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x | sigma ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = dd)
#> 
#> Distribution: gaussian1     Observations: 120
#> 
#> === mu
#> 
#> Parametric terms
#>             estimate      se     z       p  lower upper
#> (Intercept)    1.063 0.07514 14.14 < 1e-16 0.9156 1.210
#> x              1.903 0.13020 14.62 < 1e-16 1.6480 2.158
#> 
#> === sigma
#> 
#> Parametric terms
#>             estimate     se        z         p   lower   upper
#> (Intercept)  -0.9647 0.1386 -6.96100 3.375e-12 -1.2360 -0.6931
#> x             0.0216 0.2389  0.09041     0.928 -0.4466  0.4898
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -55.844138    effective df 4.00
#> cAIC 119.688    cBIC 130.838
#> fitted in 29 ms, converged
```
