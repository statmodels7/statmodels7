# The Distribution a Model Was Fitted With

The distributions7 object, with its links.

## Usage

``` r
# S3 method for class 'StatmodFit'
family(object, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

The distributions7 distribution object the model was fitted with, links
included, exactly as it was supplied. Its `@params` property names the
parameters and `@links` holds their charts.

## Details

What comes back is the family object itself, so everything the family
can do is reachable from a fit: its density, its derivatives, its
moments and its parameters' links. Compare
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), whose
[`family()`](https://rdrr.io/r/stats/family.html) returns a description.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md);
[`distributions7::distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.html)
and its siblings for what can then be asked of the family

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
family(fit)@params
#> [1] "mu"    "sigma"
```
