# The Chain Term of the Joint Gradient Along a Direction

The contributions to \\\mathrm{tr}(M\\\partial K/\partial u\[v\])\\ that
come from the recursion and not from the design: the derivative of the
term's curvature along \\v\\, the curvature re-weighted by the
derivative of the score along \\v\\, and the movement of the filter's
Jacobian, which enters through every cross derivative of the filter's
equation.

## Usage

``` r
structural_chain_extra(spec, design, jd, M, st, v)
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
  [`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md).

- v:

  The direction, over the estimated coordinates.

## Value

A single number.
