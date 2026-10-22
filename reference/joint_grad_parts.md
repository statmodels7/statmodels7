# The Direction-Free Quantities of the Joint Chain Term, Joint Filter

Computes, for a model carrying a joint filter, what
[`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md)
computes for one filter: the family's derivatives, the rows of every
equation with the joint Jacobian in place, their products with \\M\\,
the per-observation diagonals of \\M\\ and the contraction \\u\\.

## Usage

``` r
joint_grad_parts(spec, design, coef, jd, M, order = 3L)
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

- order:

  The highest order of the driving vector's derivatives read: 3 for the
  gradient, 4 for the Hessian.

## Value

A list with `u`, `V`, `Vk`, `VM`, `G`, `H`, `D3`, `D4`, `D5`, `gl`, `gw`
(the weights of the recursion, one column per filtered parameter), `bd`
(the driving data), `GF` (the positions of the filtered parameters) and
`w`.
