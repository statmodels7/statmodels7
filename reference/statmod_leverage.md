# Leverage of the Observations in a Distributional Regression

Returns, for each observation \\i\\, the block \\H\_{ii}\\ of the hat
matrix of the fit over the coefficients of all the equations. For a
model with one equation it is the usual leverage; for a model with \\d\\
equations it is a \\d \times d\\ matrix, and its trace is the number of
effective parameters that observation \\i\\ uses.

## Usage

``` r
statmod_leverage(fit, information = NULL)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- information:

  `NULL` (the default), which is
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  with the closed form of the family, or an object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).
  The default does not depend on the information that
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  resolves for the family.

## Value

A list with `blocks`, an array of dimension \\n \times d \times d\\ of
the matrices \\H\_{ii}\\, `leverage`, the \\n \times d\\ matrix of their
diagonals with the equations as columns, and `trace`, the vector of the
traces, and `clipped`, the number of observations whose information had
a negative eigenvalue set to zero (always zero with the expected
information). The \\d\\ equations are those that have coefficients.
Coordinates whose variance is missing (aliased, or held at zero by a
penalty) are not part of \\Z_i\\.

## Details

Let \\\hat\beta\\ be the stacked coefficients of the equations that have
coefficients, of length \\P\\. Write \\Z_i\\ for the \\d \times P\\
matrix whose row \\a\\ holds the row of the design of equation \\a\\ for
observation \\i\\, in the columns of that equation. The score of
observation \\i\\ in the coefficients is \\u_i = Z_i' s_i\\, where
\\s_i\\ is the vector of the derivatives of its log-likelihood in the
\\d\\ predictors (on the link scale of each parameter), multiplied by
the prior weight. Let \\W_i\\ be the \\d \times d\\ information of the
observation in those predictors, times the prior weight, and \\K\\ the
penalized information \\\sum_i Z_i' W_i Z_i + S\\ that the Bayesian
variance \\V_b = K^{-1}\\ inverts, with \\S\\ the Hessian of the
penalty. The linearized fit satisfies \\W^{-1/2}\hat s \approx (I - H)\\
W^{-1/2} s\\, with the hat matrix \$\$H = W^{1/2} Z K^{-1} Z'
W^{1/2},\$\$ where \\Z\\ stacks the \\Z_i\\, \\W\\ is block diagonal in
the \\W_i\\ and \\W^{1/2}\\ is the symmetric square root. The diagonal
block of observation \\i\\ is \$\$H\_{ii} = W_i^{1/2} Z_i K^{-1} Z_i'
W_i^{1/2}.\$\$ It is symmetric, and for an unpenalized fit with \\W\\
the information at the estimate it is a projection, so that \\\sum_i
\mathrm{tr}\\ H\_{ii} = P\\. With penalties the sum is the effective
degrees of freedom.

For a Gaussian linear model with a constant scale and the expected
information, the block of the mean equation is the scalar \\h_i =
x_i'(X'X)^{-1}x_i\\ of
[`stats::hatvalues()`](https://rdrr.io/r/stats/influence.measures.html).
For one equation of a larger model, the diagonal element \\(a, a)\\ of
the block is the leverage that
[`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md)
computes for that equation, provided that `information` is
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html),
the default. Both read \\K^{-1}\\ and not the Bayesian variance of
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which for a coordinate that the marginal criterion estimates holds the
curvature of that criterion.

The matrix \\W_i\\ and the matrix \\K\\ are both read on `information`.
With
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
the hat matrix is a projection and its trace is \\P\\ when there is no
penalty. With
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
the matrix \\W_i\\ is the observed information of the observation, which
can be indefinite; its root is then the symmetric square root of the
positive part (the eigenvalues that are not positive are set to zero,
with a zero pseudo-inverse), `clipped` counts the observations where
that happened, and the trace is close to \\P\\ but not equal to it. With
prior weights \\w_i\\ the block is the one of the row as a unit, \\w_i\\
times the block of one of its copies.

## References

Hoaglin, D. C. and Welsch, R. E. (1978). The hat matrix in regression
and ANOVA. *The American Statistician*, 32(1), 17–22.

## See also

[`hc2_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md),
which uses the blocks;
[`plot_influence()`](https://statmodels7.github.io/statmodels7/reference/plot_influence.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.3)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
lv <- statmod_leverage(fit)
range(lv$leverage[, "mu"] - hatvalues(lm(y ~ x, dd)))
#> [1] 3.469447e-18 2.498002e-16
```
