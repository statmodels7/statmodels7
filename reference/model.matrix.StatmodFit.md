# The Design of One Equation

The block of columns a distribution parameter's equation was fitted
with.

## Usage

``` r
# S3 method for class 'StatmodFit'
model.matrix(object, what = NULL, ...)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- what:

  Which distribution parameter, or `NULL` for the first.

- ...:

  Unused.

## Value

An `n x p` matrix, `p` being the width of that equation's design,
carrying the coefficient names as its column names. A Matrix object
where the equation's own blocks are sparse, a base matrix otherwise.
Zero columns for an equation carrying only a structural term.

## Details

A fit has one design per parameter, so which one is wanted is an
argument. `NULL` gives the first, as
[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
does. A term whose block moves with its coefficients is returned at the
fitted ones, which is the block the fit ended on.

A structural term contributes no columns at all, so an equation carrying
one alone gives a matrix of no columns; what such a term contributes is
a recursion, reported by
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## See also

[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for every equation at once,
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
for the coefficients these columns multiply

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + dd$x + rnorm(40, sd = 0.3)
fit <- statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
dim(model.matrix(fit))
#> [1] 40  2
colnames(model.matrix(fit, "sigma"))
#> [1] "(Intercept)" "x"          
```
