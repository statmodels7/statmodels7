# The Driving Quantity of a Filter and Its Derivatives

Replaces, in the family's link-scale derivative arrays, the entries that
a filter's recursion reads as derivatives of its driving quantity with
the derivatives of the scaled score \\u = s\\\mathcal{I}^{-d}\\.

## Usage

``` r
filter_driving(
  spec,
  theta,
  ap,
  scaling,
  gl,
  H,
  D3 = NULL,
  D4 = NULL,
  D5 = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- theta:

  The per-observation parameters.

- ap:

  The index of the parameter the filter sits in.

- scaling:

  The exponent \\d\\.

- gl, H, D3, D4, D5:

  The family's link-scale derivatives; the higher ones may be `NULL`.

## Value

A list with `gl`, `H`, `D3`, `D4` and `D5`, the same objects where
`scaling` is zero.

## Details

With \\w = \mathcal{I}^{-d}\\, Leibniz's rule over the positions of an
index tuple \\T\\ gives \$\$u_T = \sum\_{A \subseteq T} s\_{A}\\ w\_{T
\setminus A},\$\$ where \\s_A\\ is the family's derivative of order
\\1 + \|A\|\\ with the filter's parameter as one index, and Faa di
Bruno's formula over the set partitions \\\pi\\ of \\B\\ gives \$\$w_B =
\sum\_{\pi} \phi^{(\|\pi\|)}(\mathcal{I}) \prod\_{\beta \in \pi}
\mathcal{I}\_\beta, \qquad \phi^{(k)}(x) = (-d)(-d-1)\cdots(-d-k+1)\\
x^{-d-k}.\$\$ The derivatives \\\mathcal{I}\_\beta\\ are minus those of
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and its higher-order siblings, all on the link scale. The highest array
given sets the order: `D5` needs the fourth derivative of the
information.
