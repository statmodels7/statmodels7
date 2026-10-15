# The Observed Information of a Model Carrying a Mixture Term

The observed information over the coefficients of every equation
together with the own parameters of a structural term of the likelihood
shape (a regime term or a marginal break-point term), from the exact
Hessian of the mixed log-likelihood.

## Usage

``` r
statmod_regime_information(spec, ev, design, npar, offs, nb, n, w)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- ev:

  The evaluated predictors, from
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).

- design:

  The design.

- npar, offs, nb:

  The per-equation coefficient counts, their offsets and their total.

- n, w:

  The number of observations and their weights.

## Value

A square matrix over the stacked coefficients followed by the term's own
free parameters.

## Details

The matrix
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
assembles for a mixture is the complete-data information, the ordinary
one averaged over the smoothed states. It serves as a scoring matrix and
differs from the observed information: by Louis's missing-information
principle the two differ by the conditional variance of the
complete-data score, which is positive semidefinite, so the
complete-data matrix is the larger and a standard error read off it is
too small.

[`modelterms7::term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.html)
returns the exact Hessian by propagating first and second derivatives
through the same scaled forward recursion that computes the likelihood.
This function supplies the part of the calculation that belongs to the
model: how each equation's unknowns reach each predictor, and the
family's first and second derivatives at the predictor to which each
component shifts.

Unlike the filter's, these derivatives cannot be looked up from a single
evaluation: a component shifts the predictor by a level of its own, so
the family is evaluated once per component, vectorized over
observations.
