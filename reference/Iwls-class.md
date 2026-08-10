# The Iterated Weighted Least Squares Method

The S7 class of the scoring step
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
uses for the smooth block, carrying its own choice of curvature and of
decomposition.

## Usage

``` r
Iwls(
  hessian = character(0),
  approx = character(0),
  decomposition = character(0),
  maxit = integer(0),
  tol = integer(0),
  step_halving = integer(0)
)
```

## Arguments

- hessian:

  Either `"expected"` – Fisher scoring – or `"observed"`, which is
  Newton.

- approx:

  The approximation of the expected information, read only where the
  family has no closed form.

- decomposition:

  How the step is solved: `"qr"`, `"svd"`, `"chol"` or
  `"chol_crossprod"`.

- maxit:

  The iteration budget.

- tol:

  The stopping tolerance on the scaled score.

- step_halving:

  The number of halvings allowed before a step is abandoned.

## Value

An object of class `Iwls`.

## See also

[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)

## Examples

``` r
S7::S7_inherits(iwls(), Iwls)
#> [1] TRUE
```
