# Whether the Outer Derivatives Can Read a Jeffreys Penalty

Returns `TRUE` where the exact outer gradient (`order = 1`) or Hessian
(`order = 2`) can carry the movement of a Jeffreys penalty.

## Usage

``` r
jfp_outer_ok(spec, design, order = 1L)

jfp_outer_probe(d, order)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- order:

  `1` or `2`.

- d:

  A distribution.

## Value

A single logical.

## Details

The gradient reads the third derivative of the penalty, so the third
derivative of the family's expected information, and the Hessian the
fourth; each is probed once at a random point with the method that
[`jfp_emethod()`](https://statmodels7.github.io/statmodels7/reference/jfp_emethod.md)
selects. The derivatives are written for a model without a filter and
without a term whose block is recomputed at the coefficients, and the
function returns `FALSE` there.
