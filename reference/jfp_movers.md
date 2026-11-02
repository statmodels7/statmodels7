# A Recomputed Block Beside a Filter in the Jeffreys Derivatives

Returns the parts of the derivatives of the predictors that come from a
term whose block is recomputed at the coefficients
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
a smoothed
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html))
in a model with a filter: `jfp_static_rows()` the part of
\\\partial^{1+\|S\|}\eta/\partial u\\\partial S\\ for one or two
directions, one matrix per equation over every coordinate, and
`jfp_static_curv()` the part of \\\sum_t\sum_c g\_{c,t}\\
\partial^{2+\|S\|}\eta\_{c,t}/\partial u^2\partial S\\ for none, one or
two directions. `jfp_movers()` returns the terms, and `jfp_static_aug()`
the predictors' derivatives with further unknowns whose static rows are
given.

## Usage

``` r
jfp_movers(st)

jfp_static_aug(st, cols)

jfp_static_block(st, name, f)

jfp_static_pairs(st, name, f)

jfp_static_xvw(st, v, w)

jfp_static_rows(st, dirs)

jfp_static_curv(st, g, dirs = list())
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md),
  of a model with a filter.

- cols:

  A list of column blocks, each a list with `name` and `rows` (one
  matrix or `NULL` per distribution parameter).

- name:

  The name of a block of columns.

- f:

  A function of a unit returning the block's static rows in the unit's
  equation.

- v, w:

  Directions over every coordinate.

- dirs:

  A list of directions over every coordinate (the full width of the
  filter's rows).

- g:

  A named list of per-observation weights, already multiplied by the
  observation weights.

## Value

`jfp_static_rows()` one matrix per distribution parameter over every
coordinate; `jfp_static_curv()` a square matrix over every coordinate;
`jfp_movers()` the units of
[`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md);
`jfp_static_aug()` what
[`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md)
returns, with `pos` (the columns of each block), `oi` (the columns of
the original coordinates), `pad()` (a direction over the original
coordinates placed on the wider width) and `unit()` (the direction of
one added column).

## Details

A filter writes \\\eta = F(s, \psi)\\ from the static predictors
\\s(u)\\ of every equation and its own parameters. The kernels of the
filter
([`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md))
differentiate \\F\\ with the static rows \\R = \partial s/\partial u\\
held at their value, which is the whole derivative only where \\s\\ is
linear in \\u\\. Where a block moves, the derivatives of \\s\\ of order
two and more are \\s''\[a, b\] = \partial R\[a\]\\b\\ and so on
([`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)
and its siblings), and the derivative of \\\eta\\ along a set of indices
is the sum over the partitions of the set into blocks, a block of one
index entering \\F\\ as a full direction and a block of \\j \ge 2\\
indices as the static perturbation \\s^{(j)}\\ along them. The kernels
return every such term once the perturbations are placed as further
unknowns: a column with the static rows \\x\\ and no part in \\\psi\\ is
a direction along which \\F\\ moves by \\x\\ alone, so its entry in the
kernels' result is the derivative of \\F\\ with \\x\\ in one argument. A
free index that lies in a block of two or more indices is carried by one
such column per value of the index.

The partitions written out are those of \\\\k, v\\\\ and \\\\k, v, w\\\\
for the rows, and of \\\\k, l\\\\, \\\\k, l, v\\\\ and \\\\k, l, v,
w\\\\ for the contraction, without the partition into singletons, which
is the kernels' own result. In an equation that no filter drives \\F\\
is the identity in that equation's static predictor, so the unknown's
column enters there only through the contraction of the Jacobian.
