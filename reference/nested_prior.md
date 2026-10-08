# The Prior of One Nested Random-Effect Term

The Gaussian prior of a random-effect sub-term, as
[`random_prior()`](https://statmodels7.github.io/statmodels7/reference/random_prior.md)
returns it for a term written in an equation.

## Usage

``` r
nested_prior(spec, design, fit, row)
```

## Arguments

- spec, design:

  The specification and its design.

- fit:

  The fitted model, for the hyperparameters.

- row:

  One row of
  [`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md).

## Value

A list with `gaussian = TRUE`, `dim` and `chol`.

## Details

The penalty unit is the outer term's, keyed by the outer term where it
carries one penalty and by the row's key where it carries several. Only
a Gaussian prior is integrated: inside a nonlinear term the average is
taken over the nodes of a product grid at which the term is evaluated,
and a prior of another family signals an error naming `random = "zero"`.
