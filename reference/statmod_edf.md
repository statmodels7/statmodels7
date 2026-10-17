# Effective Degrees of Freedom, Per Term

Computes the effective degrees of freedom of each term: from the
diagonal of the model's smoother where a penalized term is not kinked,
from the number of coefficients away from the kink for a kinked penalty,
from the term's free parameters for a structural term, and through
modelterms7's `edf()` for the remaining terms.

## Usage

``` r
statmod_edf(
  spec,
  coef,
  design,
  hyper,
  expected = TRUE,
  approx = "opg",
  aliased = integer(0)
)
```

## Arguments

- spec:

  The specification.

- coef:

  The coefficients.

- design:

  The design.

- hyper:

  The hyperparameters.

- expected:

  Whether the curvature is the expected information.

- approx:

  The approximation for the expected information.

- aliased:

  Stacked positions of the coefficients that the fit did not estimate;
  they are not counted.

## Value

A data frame with one row per term and the columns `parameter`, `term`,
`coefficients` and `edf`, or `NULL` where the model has no term. `edf`
is `NA` for a term whose count could not be read.

## Details

A smooth penalized term takes its share from the diagonal of
\\(H+S)^{-1}H\\ for the whole model, read over the coefficients of every
equation at once. The computation needs the unpenalized curvature of
every equation, and not only the coefficients and hyperparameters of the
term. That curvature is computed once for all terms and not per term.

A term may carry more than one penalty, over different parameters of its
own, and each has hyperparameters of its own filed under a key of its
own. The result has one row per term, the granularity of a table of
terms, and the hyperparameters are handed over keyed by the penalty
names that
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
returns, which is the shape in which `edf()` reads them. Passing the
hyperparameters of one penalty for a term carrying two would count the
whole block against it.

The arguments of `edf()` are passed by name: its third argument is the
curvature and its fourth the hyperparameters.
