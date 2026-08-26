# A Refreshable Block's Second Derivative Along Two Directions

[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html)
on one unit, with the shape it returns checked rather than assumed.

## Usage

``` r
refresh_dblock2(un, v, u, n)
```

## Arguments

- un:

  One entry of
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- v, u:

  The two directions, each as long as the term's coefficients.

- n:

  The number of observations.

## Value

A numeric matrix, `n` by the term's coefficient count.

## Details

A term that has not written the method inherits the base one and gets
zeros, which is exact for a design that does not move and a deliberate
refusal for a block that is a working linearization rather than a
Jacobian. Nothing is differenced here in either case.

## See also

[`refresh_dblock()`](https://statmodels7.github.io/statmodels7/reference/refresh_dblock.md),
[`trace_refresh4()`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md)
