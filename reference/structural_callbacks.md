# The Score and Curvature a Filter Is Driven By

The derivative of the log-density in one distribution parameter's
predictor, and its second derivative, as functions of that predictor at
one observation, with every other parameter held where it is.

## Usage

``` r
structural_callbacks(spec, theta, p)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- theta:

  The per-observation parameters, as
  [`statmod_eta`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- p:

  The distribution parameter the term sits in.

## Value

A list with `score`, `curvature` and `logdens`.

## Details

They are read one observation at a time because a score-driven filter
evaluates them at the predictor it has just produced, which is not known
before the recursion reaches that observation. A term whose callbacks do
not depend on the state – a regime chain, whose levels shift a predictor
known in advance – is given the whole index vector instead.
