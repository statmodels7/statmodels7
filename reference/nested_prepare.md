# A Prediction With Nested Random Effects Set Aside

Prepares the specification, the design and the coefficients of a
prediction in which random effects written inside a subformula are read
at zero or averaged over their prior.

## Usage

``` r
nested_prepare(fit, spec, nest, unseen)
```

## Arguments

- fit:

  The fitted model.

- spec:

  The specification at the prediction's rows.

- nest:

  The rows of
  [`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md)
  set aside, with `mode`.

- unseen:

  As in
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md).

## Value

A list with `spec`, `design`, `coef` and `cols`, the latter from
[`nested_columns()`](https://statmodels7.github.io/statmodels7/reference/nested_columns.md).

## Details

Every term holding such an effect is replaced by the copy
[`pool_nested_random()`](https://statmodels7.github.io/statmodels7/reference/pool_nested_random.md)
returns, its effects' coefficients are set to zero and its columns are
recomputed there. The columns of the effects themselves are then set to
zero, so that the delta method of the prediction carries the uncertainty
of the other coefficients alone, as it does for a term written in an
equation.
