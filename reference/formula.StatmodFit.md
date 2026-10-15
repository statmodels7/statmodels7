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
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- ...:

  Unused.

## Value

The formula object supplied to
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
with its bars and its original environment intact.

## Details

The formula is returned whole, bars included, because the bars are part
of what was written.
[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
splits the formula into one equation per distribution parameter; it
takes the formula and the names of the parameters.

## See also

[`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
to split it into one equation per parameter,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
formula(fit)
#> y ~ x | sigma ~ x
#> <environment: 0x5603fe82c298>
```
