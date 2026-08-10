# Summarize a Fitted Model

One coefficient table per distribution parameter – estimate, standard
error, Wald statistic, p-value and interval – with the degrees of
freedom, the information criteria and the qualifications the numbers
carry.

## Usage

``` r
# S3 method for class 'StatmodFit'
summary(object, level = 0.95, type = c("bayesian", "frequentist"), ...)
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

**A hyperparameter carries no standard error yet.** It is held at the
value it was given rather than estimated, so the row reports the value
and marks it fixed; inventing an interval for a number nothing estimated
would be worse than the empty column. Estimating them by an outer
criterion is what fills those rows in.

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
#> (Intercept) -0.96470 0.1386 -6.96100 3.375e-12 -1.2360 -0.6931
#> x            0.02161 0.2389  0.09044    0.9279 -0.4466  0.4898
#> 
#> 95% intervals, bayesian variance
#> log-likelihood -55.844138    df 4.00    AIC 119.688    BIC 130.838
#> fitted in 31 ms, converged
```
