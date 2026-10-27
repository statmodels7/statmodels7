# The Working Response and Weights of One Equation

Computes the working weights \\h_i\\ and the working response \\z_i =
\eta_i + s_i/h_i\\, which define the weighted least squares problem that
locally represents the log-likelihood.

## Usage

``` r
coord_working(spec, ep, coef, design, p, expected, emethod)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- ep:

  The linear predictors and the parameters that they imply, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- coef:

  A named list of coefficient vectors.

- design:

  The design.

- p:

  Which distribution parameter's equation, a string.

- expected:

  `TRUE` for the expected information, `FALSE` for the observed one.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

A list with `w` and `z`, each a numeric vector of length `spec@n_obs`. A
row of weight zero has `w` zero and `z` equal to its predictor. The
value is `NULL` where the working problem is not usable: a non-finite or
non-positive \\h_i\\, or a non-finite score, on a row of positive
weight. The observed information can produce non-positive values far
from the optimum.

## Details

For a Gaussian response on the identity link the quadratic is exact and
the first iteration solves the problem. Elsewhere it is the local
approximation on which a scoring step works, and the weights are rebuilt
at each iteration. A row of weight zero has no curvature and no score,
and it enters the descent with weight zero at its own predictor.

Where the design carries a structural term, the predictor of an equation
includes a part that belongs to that term, and the columns of a block do
not carry it. The working response is then built on the static
predictor, with \\s_i\\ the score of the model from
[`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md),
so the quadratic has the gradient of the log-likelihood at the current
coefficients. The weights stay the family's curvature at the whole
predictor. Another positive weight would change the path of the descent
and would leave its fixed point where it is.
