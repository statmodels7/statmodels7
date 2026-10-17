# The Covariance of Every Predictor at One Row, for a New Group

Returns, for each observation, the covariance over the distribution's
parameters of their linear predictors: the estimation uncertainty of the
fixed part, from the variance matrix, plus the variance that a new
group's effects add.

## Usage

``` r
predictive_cov(object, spec, design, blocks, fixed = TRUE, ...)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  The specification and design at the rows predicted.

- blocks:

  The blocks of
  [`random_blocks()`](https://statmodels7.github.io/statmodels7/reference/random_blocks.md).

- fixed:

  `FALSE` to skip the fixed part, which is then zero.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

A list with `fixed` and `random`, arrays of dimension `P x P x n` over
the parameters.

## Details

The fixed part is \\x\_{ip}^\top V\_{pq} x\_{iq}\\, evaluated after the
columns of the effects set aside are set to zero, so it is the
uncertainty of the typical group's predictor; the part that the effects
add is \\z\_{ip}^\top \Sigma\_{pq} z\_{iq}\\, summed over the blocks of
[`random_blocks()`](https://statmodels7.github.io/statmodels7/reference/random_blocks.md).
