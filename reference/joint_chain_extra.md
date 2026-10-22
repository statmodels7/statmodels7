# The Chain Term of the Joint Gradient Along a Direction, Joint Filter

Computes, for a model carrying a joint filter, the three contributions
that
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
computes for one filter: the third derivative of the levels along \\v\\
contracted against the weighted score, the curvature of the levels
re-weighted by the derivative of the score along \\v\\, and the movement
of the joint Jacobian through every cross derivative of a filtered
equation.

## Usage

``` r
joint_chain_extra(spec, design, jd, M, st, v)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows.

- M:

  The matrix against which the trace is taken.

- st:

  The shared quantities, from
  [`joint_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/joint_grad_parts.md).

- v:

  The direction, over the estimated coordinates.

## Value

A single number.
