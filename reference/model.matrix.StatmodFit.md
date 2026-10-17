# The Design of One Equation

Returns the block of columns used to fit the equation of a distribution
parameter.

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

An `n x p` matrix, where `p` is the width of the design of that
equation, with the coefficient names as its column names. The result is
a Matrix object where the blocks of the equation are sparse and a base
matrix otherwise. The matrix has zero columns for an equation that
carries only a structural term.

## Details

A fit has one design per parameter, and `what` selects the parameter.
`NULL` selects the first, as
[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
does. A term whose block moves with its coefficients is returned at the
fitted coefficients, which is the block on which the fit ended.

A structural term contributes no columns at all, so an equation that
carries only a structural term gives a matrix of no columns. The
contribution of such a term is a recursion, reported by
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md).

## See also

[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for every equation at once,
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
for the coefficients that these columns multiply

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
