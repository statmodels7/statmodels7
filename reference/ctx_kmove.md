# How the Penalized Matrix Moves With the Coefficients

The array
[`u_vector`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
contracts against the leverage diagonal, together with the builder that
reads it: the third derivative of the log-density where the criterion
uses the observed information, and the derivative of the expected
information where it uses the expected one.

## Usage

``` r
ctx_kmove(ctx, spec, design, coef, hyper, method)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

A list with `deriv` and `key`.

## Details

\\K\\ enters the criterion through its determinant, so the gradient
needs \\\partial K/\partial\beta\\. With the observed information that
is \\-\ell'''\\, which every family carries in closed form. With the
expected one it is \\-\partial\\\mathbb{E}\[\ell''\]/\partial\eta\\,
which is NOT \\-\mathbb{E}\[\ell'''\]\\: differentiating an expectation
moves the measure as well as the integrand, and the missing piece
\\\mathbb{E}\[\ell\_{ab}\ell\_{c}\]\\ is a mixed moment no Bartlett
identity isolates – the third ties the symmetrized sum, not the single
term.

distributions7 supplies it as
[`distrib_dexpected_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html).
The two arrays are keyed differently, the observed one being symmetric
in all three indices and the expected one in its first two only, so the
key builder travels with the array rather than being assumed by the
consumer.

## See also

[`u_vector`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`outer_gradient_ok`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
