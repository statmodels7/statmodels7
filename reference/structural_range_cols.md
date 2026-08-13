# Which of a Structural Term's Free Parameters a Penalty Covers

Positions among the term's free parameters that some penalty of the
model shrinks, which are the directions
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
integrates over.

## Usage

``` r
structural_range_cols(spec, design, key, free)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- key:

  The term's name.

- free:

  The term's free parameters, in order.

## Value

An integer vector.
