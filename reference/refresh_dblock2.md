# A Refreshable Block's Second Derivative Along Two Directions

Calls
[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html)
on one unit and checks the shape of the result. A block of another shape
signals an error that names the term.

## Usage

``` r
refresh_dblock2(un, v, u, n)
```

## Arguments

- un:

  One entry of
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- v, u:

  The two directions, each of which is as long as the term's
  coefficients.

- n:

  The number of observations.

## Value

A numeric matrix, `n` by the term's coefficient count.

## Details

The base method on `model_term` returns zeros, which is exact for a
design that does not move. A sharp
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
returns zeros because its second derivative vanishes away from the
break-points, and a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
returns zeros because its block is a working linearization with a frozen
weight and not a Jacobian. The function does not difference in any case.

## See also

[`refresh_dblock()`](https://statmodels7.github.io/statmodels7/reference/refresh_dblock.md),
[`trace_refresh4()`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md)
