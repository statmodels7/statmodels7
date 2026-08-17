# How the Determinant Reads a Block That Moves With the Coefficients

The part of \\u_c = \mathrm{tr}(M\\\partial K/\partial\beta_c)\\ that
[`u_vector`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
does not compute: everything coming from \\\partial X/\partial\beta\\
where a term's block depends on its own coefficients.

## Usage

``` r
u_refresh(
  spec,
  design,
  coef,
  M,
  params,
  npar,
  offs,
  total,
  expected = FALSE,
  approx = "bartlett",
  units = NULL,
  Hl = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, already refreshed at `coef`.

- coef:

  The coefficients at the penalized mode.

- M:

  The matrix the trace is taken against.

- params, npar, offs, total:

  The block bookkeeping.

- expected:

  Whether the criterion carries the expected information.

- approx:

  The approximation for the expected information.

- units:

  The refreshable terms, from
  [`refresh_units`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md),
  or `NULL` to resolve them here.

- Hl:

  The link-scale curvature, from
  [`refresh_hessian`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md),
  or `NULL` to compute it here.

## Value

A numeric vector as long as the stacked coefficients.

## Details

With \\H\_{(a,j),(b,k)} = -\sum_i w_i \ell\_{ab,i}X_a\[i,j\]X_b\[i,k\]\\
and a block that moves, differentiating in \\\beta_c\\ gives three
contributions and
[`u_vector`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
computes one. The other two are transposes under the trace, so with
\\R\_{ab}\[i,j\] = \sum_k M\_{(a,j),(b,k)}X_b\[i,k\]\\ and \\A_a\[i,j\]
= w_i\sum_b \ell\_{ab,i}R\_{ab}\[i,j\]\\, \$\$\Delta u_c =
-2\sum\_{i,j}A_a\[i,j\]\\\partial X_a\[i,j\]/\partial\beta_c.\$\$

**The derivative is asked of the TERM**, through
[`term_block_contract`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html),
and never differenced here. Two reasons, both measured: a term knows its
own chain rule – the links on its parameters and a subformula's design –
and a break-point column is a step function in its break-point, so a
difference quotient of it diverges as the step shrinks rather than
converging. A term that does not implement the contraction inherits
zeros, which is exactly right for a fixed design.

## See also

[`u_vector`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`refresh_amat`](https://statmodels7.github.io/statmodels7/reference/refresh_amat.md),
[`term_block_contract`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.html)
