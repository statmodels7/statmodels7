# The Penalty of a Specification at Given Coefficients

The sum of every penalized term's `penalty_value()`, with its gradient
and Hessian placed in the columns that term owns.

## Usage

``` r
statmod_penalty_at(
  spec,
  coef,
  hyper,
  design = statmod_design(spec),
  what = c("value", "gradient", "hessian")
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- hyper:

  A named list, one entry per parameter, each a named list of
  hyperparameter vectors per penalized term.

- design:

  The design.

- what:

  One of `"value"`, `"gradient"`, `"hessian"`.

## Value

A number, a named list of vectors, or a square matrix.
