# Simulate a Response From a Written Model

Takes a formula, a distribution and a data frame of covariates, and
draws a response from that model – at coefficients the caller supplies,
or at coefficients drawn at random.

## Usage

``` r
rstatmod(formula, distrib, data, par = NULL, sd = 1, offsets = NULL)
```

## Arguments

- formula:

  The model formula, as
  [`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes it.

- distrib:

  A distributions7 distribution object.

- data:

  A data frame of covariates.

- par:

  Optional named list of coefficient vectors.

- sd:

  The standard deviation of the drawn coefficients.

- offsets:

  Optional named list of offsets.

## Value

The data frame with the response added, carrying the attributes `"par"`
(the coefficients used) and `"theta"` (the parameters they gave).

## Details

This is not [`simulate`](https://rdrr.io/r/stats/simulate.html), which
draws from a model that has already been fitted. The `r` prefix is R's
own for a random draw, so the two names cannot be confused.

The point of it is to have data whose truth is known: write the model,
draw from it, fit it back, and see whether the fit recovers what was put
in. A covariate needs no declaring – a factor becomes its contrasts and
a numeric stays itself, because the design comes from the same
interpreter a fit uses.

**The coefficients.** `par = NULL` draws them, each independently from
`rnorm(1, 0, sd)`, which on the link scale gives predictors of order
one. A named list fixes them instead, one vector per distribution
parameter in the design's order; a parameter left out of that list is
drawn. [`coef()`](https://rdrr.io/r/stats/coef.html) on the result
reports what was used, drawn or given, so a simulation is reproducible
from its own output.

**The response's name** is the formula's left-hand side when it is a
symbol, and `"y"` otherwise.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))

# coefficients drawn
sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
attr(sim, "par")
#> $mu
#> (Intercept)           x          gb 
#> -0.05612874 -0.15579551 -1.47075238 
#> 
#> $sigma
#> (Intercept) 
#>  -0.4781501 
#> 

# or given, and recovered by a fit
sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                 par = list(mu = c(1, 2), sigma = log(0.3)))
statmod(y ~ x, distributions7::gaussian1_distrib(), sim2)
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ x, distrib = distributions7::gaussian1_distrib(), 
#>             data = sim2)
#> 
#> Distribution: gaussian1
#> Observations: 50
#> 
#>   mu         ~ x
#>                linpar           2 coef
#>   sigma      ~ 1
#>                linpar           1 coef
#> 
#> log-likelihood 0.466254    objective -0.466254
#> fitted in 29 ms, converged
```
