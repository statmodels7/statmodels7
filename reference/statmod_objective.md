# The Objective, Its Gradient and Its Hessian, Stacked

\\F(\beta) = -\ell(\beta) + \sum_t \rho_t\\, unaveraged, over the
coefficients of every parameter stacked into one vector.

## Usage

``` r
statmod_objective(
  spec,
  hyper,
  design = statmod_design(spec),
  expected = TRUE,
  approx = "bartlett"
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- hyper:

  The hyperparameters.

- design:

  The design.

- expected:

  Whether the information is the expected one.

- approx:

  The approximation for the expected information.

## Value

A list of functions `fn`, `gr` and `he` of the stacked coefficient
vector, plus `split` and `stack` to move between that vector and the
per-parameter list.

## Details

The objective is NOT divided by the sample size. A penalty is a negative
log-prior at full size, and a posterior is a log-likelihood plus a
log-prior at full size; averaging the first would make a hyperparameter
mean something that depends on \\n\\. What is scaled instead is the
stopping rule, in the one place it is read.

## See also

[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
