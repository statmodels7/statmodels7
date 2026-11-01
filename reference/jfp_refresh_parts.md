# The Movement of a Jeffreys Block Inside a Refreshed Term

Returns, for a block inside a subformula of a term whose block is
recomputed at the coefficients
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)),
the term at the coefficients and the derivative of each of the block's
basis columns in the term's coefficients.

## Usage

``` r
jfp_refresh_parts(st, b)
```

## Arguments

- st:

  The state under construction, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- b:

  A
  [`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
  of kind `"refresh"`.

## Value

A list with `tm`, `bt`, `tcols` (the term's columns in its equation),
`tidx` (their positions in the joint vector), `pos` (the basis columns
among the term's) and `D` (one matrix per basis column, observations by
the term's coefficients).

## Details

The block's columns are columns of the term's Jacobian, so \\\partial
X\_{t,i}/\partial\beta_l = \partial^2\eta_t/\partial
\beta_i\partial\beta_l\\, which is column \\l\\ of the block's
derivative along \\e_i\\,
[`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html).
One call per basis column gives every derivative the gradient and the
Hessian read.
