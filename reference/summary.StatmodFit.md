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

**What a Wald p-value means here depends on the row**, and the summary
says which is which rather than printing one column and leaving it at
that. For an unpenalized coefficient it is the usual thing. For a
coefficient in a penalized block it is conditional on the smoothing
parameter, which was not estimated jointly with it, and it does not
account for the shrinkage of the estimate towards zero. For a block a
kinked penalty selected – a lasso, a SCAD, an MCP – the row exists only
because that coefficient survived the selection, and a naive interval
there under-covers; the coefficients set exactly to zero carry `NA`,
since at the kink there is no curvature to read.

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
#> mu
#>                      estimate      se     z       p  lower upper  
#> linpar / (Intercept)    1.063 0.07514 14.14 < 1e-16 0.9156 1.210  
#> linpar / x              1.903 0.13020 14.62 < 1e-16 1.6480 2.158  
#> 
#> sigma
#>                      estimate     se        z         p   lower   upper  
#> linpar / (Intercept) -0.96470 0.1386 -6.96100 3.375e-12 -1.2360 -0.6931  
#> linpar / x            0.02161 0.2389  0.09044    0.9279 -0.4466  0.4898  
#> 
#> 95% intervals, bayesian variance
#> 
#> log-likelihood -55.844138    df 4.00    AIC 119.688    BIC 130.838
#> fitted in 29 ms, converged
```
