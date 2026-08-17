# A Refreshable Block's Derivative Along One Direction

[`term_block_deriv`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)
on one unit, with the shape it returns checked rather than assumed.

## Usage

``` r
refresh_dblock(un, v, n)
```

## Arguments

- un:

  One entry of
  [`refresh_units`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- v:

  The direction, as long as the term's coefficients.

- n:

  The number of observations.

## Value

A numeric matrix, `n` by the term's coefficient count.
