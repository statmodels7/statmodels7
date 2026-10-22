# The Derivative of the Joint Matrix Along a Direction, Joint Filter

Computes, for a model carrying a joint filter, the matrix that
[`structural_dk_matrix()`](https://statmodels7.github.io/statmodels7/reference/structural_dk_matrix.md)
computes for one filter: the part of \\\partial K/\partial u\[v\]\\ that
comes from the family and the recursion.

## Usage

``` r
joint_dk_matrix(spec, design, jd, st, v)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows.

- st:

  The shared quantities, from
  [`joint_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/joint_grad_parts.md).

- v:

  The direction, over the estimated coordinates.

## Value

A square symmetric matrix over the estimated coordinates.

## Details

Each term of the single-filter matrix that reads the filter's equation
is a sum over the filtered parameters here. The third derivative of the
levels and their re-weighted curvature come from one call each to the
joint recursion, with one column of weights per filtered parameter, and
the movement of every filtered row enters through the cross derivatives
of its own parameter.
