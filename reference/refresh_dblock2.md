# A Refreshable Block's Second Derivative Along Two Directions

[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html)
on one unit, with the shape of the result checked. A block of another
shape signals an error that names the term.

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

The base method on `model_term` answers zeros, which is exact for a
design that does not move. A sharp
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
answers zeros because its second derivative vanishes away from the
break-points, and a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
answers zeros because its block is a working linearization with a frozen
weight and not a Jacobian. Nothing is differenced here in any case.

## See also

[`refresh_dblock()`](https://statmodels7.github.io/statmodels7/reference/refresh_dblock.md),
[`trace_refresh4()`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md)
