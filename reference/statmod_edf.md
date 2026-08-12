# Effective Degrees of Freedom, Per Term

Asks each term what it spends, through modelterms7's `edf()`.

## Usage

``` r
statmod_edf(spec, coef, design, hyper, expected = TRUE, approx = "bartlett")
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

## Value

A data frame with one row per term.

## Details

A smooth penalized term counts \\\mathrm{tr}\[(H+S)^{-1}H\]\\ over its
own block, so it needs the unpenalized curvature there and not only its
coefficients and its hyperparameters. That block is cut out of the
likelihood's information, which is computed once for every term rather
than per term.

A term may carry more than one penalty, over different parameters of its
own, and each has hyperparameters of its own filed under a key of its
own. The row stays per term, which is the granularity a table of terms
wants, and the hyperparameters are handed over keyed by the penalty
names
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html)
gives, which is the shape `edf()` reads them in. Passing the
hyperparameters of one penalty for a term carrying two would count the
whole block against it.

The arguments are passed BY NAME. `edf()`'s third argument is the
curvature and its fourth the hyperparameters, and a positional call put
the hyperparameters where the curvature belongs: every smooth term then
reported `NA`, the total degrees of freedom counted the unpenalized
terms alone, and AIC and BIC were built on that count.
