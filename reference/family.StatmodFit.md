# The Distribution of a Fitted Model

Returns the distributions7 object of the model, with its links.

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

The distributions7 distribution object used to fit the model, including
its links. A binomial family fitted to a two-column response carries the
row sums as its `size`, so in that case the object differs from the one
supplied. Its `@params` property names the parameters and `@link_params`
holds their links.

## Details

The object returned is the family itself, so its density, its
derivatives, its moments and its parameters' links are available from a
fit. In contrast, [`stats::glm()`](https://rdrr.io/r/stats/glm.html) has
a [`family()`](https://rdrr.io/r/stats/family.html) that returns a
description.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md);
[`distributions7::distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.html)
and its siblings for the quantities that can be computed from the family

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
family(fit)@params
#> [1] "mu"    "sigma"
```
