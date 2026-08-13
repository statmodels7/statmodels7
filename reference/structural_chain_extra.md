# The Two Pieces of the Chain Term That Read the Direction

The contributions to \\\mathrm{tr}(M\\\partial K/\partial u\[v\])\\ that
come from the recursion rather than from the design: the derivative of
the filter's own Jacobian, and the derivative of the term the level
contributes to the information.

## Usage

``` r
structural_chain_extra(spec, design, jd, M, st, v)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows.

- M:

  The matrix the trace is taken against.

- st:

  The shared quantities, from
  [`structural_grad_parts`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md).

- v:

  The direction, over the estimated coordinates.

## Value

A single number.
