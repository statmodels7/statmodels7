# The Score and Curvature a Filter Is Driven By

The derivative of the log-density in one distribution parameter's
predictor, and its second derivative, as functions of that predictor at
one observation, with every other parameter held where it is.

## Usage

``` r
structural_callbacks(spec, theta, p, scaling = 0)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- theta:

  The per-observation parameters, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- p:

  The distribution parameter the term sits in.

- scaling:

  The exponent \\d\\ of
  [`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html):
  where it is not zero, `score` returns \\u = s\\\mathcal{I}^{-d}\\ and
  `curvature` its derivative in the predictor, \\u' =
  s'\\\mathcal{I}^{-d} - d\\s\\\mathcal{I}^{-d-1}\mathcal{I}'\\, and the
  compiled context carries `scaling`, its kernel composing the same
  quantities from the family's scalar entry points.

## Value

A list with `score`, `curvature` and `logdens`.

## Details

They are read one observation at a time because a score-driven filter
evaluates them at the predictor it has just produced, which is not known
before the recursion reaches that observation. A term whose callbacks do
not depend on the state is given the whole index vector instead. A
regime chain is such a term: its levels shift a predictor known in
advance.
