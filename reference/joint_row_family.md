# The Family at One Observation, in Every Predictor

Returns a function of the row and of the predictors of every equation at
that row, which gives the link-scale gradient and Hessian of the
log-density there and, where `info` is `TRUE`, the expected information
and its derivative in each predictor.

## Usage

``` r
joint_row_family(spec, info = FALSE)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- info:

  Whether the expected information and its derivative are read.

## Value

A function `function(i, e, yi)` of the row, of the vector of all the
predictors at that row and of the response there (by default the
specification's own), returning a list with `g`, `H` and, where `info`
is `TRUE`, `I` and `dI` (an array whose slice `[, , c]` is the
derivative of `I` in the predictor of parameter `c`).

## Details

A joint filter reads its driving vector at the predictor that it has
just produced, one observation at a time, so the family is evaluated at
one point per call. The S7 methods of the family and of its links, and
the names of the components, are resolved once, and the chain rule onto
the link scale is applied here: with \\h_a\\ the inverse link of
parameter \\a\\, the gradient is \\\ell_a h_a'\\ and the Hessian
\\\ell\_{ab} h_a' h_b' + \delta\_{ab}\\\ell_a h_a''\\; the expected
information \\\mathcal{I}\_{ab} = -\mathbb{E}\[\ell\_{ab}\]\\h_a' h_b'\\
carries no \\h''\\ term, because the score has mean zero, and its
derivative in the predictor of parameter \\c\\ is \$\$-\frac{\partial
\mathbb{E}\[\ell\_{ab}\]}{\partial\theta_c}\\ h_c' h_a' h_b' -
\mathbb{E}\[\ell\_{ab}\]\\(\delta\_{ac} h_a'' h_b' + \delta\_{bc} h_a'
h_b''),\$\$ as
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
composes it. Where the family's expected information is not written out
([`distributions7::expected_hessian_exact()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_exact.html)
is `FALSE`) or its derivative is the default difference, the two
generics are called on the link scale instead, so that their own
approximation is the one read.
