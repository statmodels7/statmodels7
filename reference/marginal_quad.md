# A Marginal Prediction by Adaptive Quadrature Over One Univariate Prior

The average over a new group's effect when the only effect set aside is
one coordinate with a prior that is not Gaussian: each quantity is
\\\int q(\eta_0 + z u)\\ f_b(u)\\ du\\, with \\f_b\\ the prior's own
density, computed by
[`numericals7::quad_vec()`](https://statmodels7.github.io/numericals7/reference/quad_vec.html)
over the rows at once.

## Usage

``` r
marginal_quad(spec, eta0, Z, param, prior, what, kind, inv)
```

## Arguments

- spec:

  The specification at the rows predicted.

- eta0:

  The predictors at the fixed part, a named list.

- Z:

  The within-group column of the term, one value a row.

- param:

  The parameter whose equation holds the term.

- prior:

  What
  [`random_prior()`](https://statmodels7.github.io/statmodels7/reference/random_prior.md)
  returns for the term.

- what, kind:

  What was asked for and its kind, as in
  [`random_marginal()`](https://statmodels7.github.io/statmodels7/reference/random_marginal.md).

- inv:

  The inverse link or its derivative.

## Value

As
[`random_marginal()`](https://statmodels7.github.io/statmodels7/reference/random_marginal.md).

## Details

Measured on a Bernoulli with a logit link and a Student t prior, 10000
Monte Carlo draws were 3e-3 out against
[`integrate()`](https://rdrr.io/r/stats/integrate.html); this is within
the default tolerances of `quad_vec()`.
