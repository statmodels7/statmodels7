# What the Joint Chain Term Needs Before a Direction Is Known

The quantities of
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md)
that do not depend on which hyperparameter is being differentiated: the
family's derivatives, the filter's forward Jacobian, the per-observation
diagonal of \\M\\ and the contraction \\u\\.

## Usage

``` r
structural_grad_parts(spec, design, coef, jd, M)

structural_grad_parts_impl(spec, design, coef, jd, M)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- jd:

  The joint rows, from
  [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md).

- M:

  The matrix against which the trace is taken.

## Value

A list with `u`, `V`, `Vk`, `VM`, `G`, `H`, `D3`, `D4`, `gl`, `s_at`,
`sd_at`, `cd_at`, `seed`, `blocks`, `blocks_data` and `w`: the
contraction \\u\\, the rows of the joint predictor derivative with the
filter's Jacobian in place, their products with \\M\\, the
per-observation diagonals of \\M\\, the family's derivatives, the score
and curvature that drive the recursion, the seed of the recursion, the
`blocks` callback with its data, and the observation weights.

## Details

Computing them once is what keeps a model with several hyperparameters
affordable: only the two pieces that read the direction are repeated,
and each of those is one pass of the recursion.
