# The Fitted Distribution of the Effects of a Random Term

Reads the penalty of the term and the fitted hyperparameters. Returns
the marginal density of each coefficient within a level and, for a
gaussian prior, the fitted covariance matrix of the coefficients.

## Usage

``` r
random_prior_fit(fit, param, name, tm)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param, name:

  The equation and the name of the term.

- tm:

  The random term.

## Value

A list with `wnames`, `dens` (a list of functions, one for each
coefficient within a level, `NULL` where the density is not available),
`sd` (the standard deviation of the gaussian prior of each coefficient,
`NA` otherwise) and `Sigma` (the covariance matrix, `NULL` unless the
prior is gaussian).
