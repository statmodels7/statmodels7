# Whether the Outer Hessian Reads Every Recomputed Block

Returns `TRUE` where every term that recomputes its block is one for
which
[`structural_refresh_second()`](https://statmodels7.github.io/statmodels7/reference/structural_refresh_second.md)
has been verified:
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
smoothed or not. Other terms answer `FALSE`.

## Usage

``` r
refresh_order2_ok(design)
```

## Arguments

- design:

  The design.

## Value

A single logical.
