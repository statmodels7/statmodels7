# The Formula a Model Was Written With

The formula as supplied, every distribution parameter's equation
included.

## Usage

``` r
# S3 method for class 'StatmodFit'
formula(x, ...)
```

## Arguments

- x:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

A formula.

## Details

It is returned whole rather than split into one formula per parameter:
the bars are part of what was written, and a caller wanting the
equations separately gets them from the fit's specification.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
formula(fit)
#> y ~ x | sigma ~ x
#> <environment: 0x55e593daf8a0>
```
