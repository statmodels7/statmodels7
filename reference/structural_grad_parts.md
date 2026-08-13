# What the Joint Chain Term Needs Before a Direction Is Known

The quantities of
[`statmod_structural_grad`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md)
that do not depend on which hyperparameter is being differentiated: the
family's derivatives, the filter's forward Jacobian, the per-observation
diagonal of \\M\\ and the contraction \\u\\.

## Usage

``` r
structural_grad_parts(spec, design, coef, jd, M)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- jd:

  The joint rows, from
  [`joint_design_rows`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md).

- M:

  The matrix the trace is taken against.

## Value

A list of the shared quantities.

## Details

Computing them once is what keeps a model with several hyperparameters
affordable: only the two pieces that read the direction are repeated,
and each of those is one pass of the recursion.
