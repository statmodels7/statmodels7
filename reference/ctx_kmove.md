# How the Penalized Matrix Moves With the Coefficients

Returns the array that
[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)
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
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

A list with `deriv` and `key`.

## Details

\\K\\ enters the criterion through its determinant, so the gradient
needs \\\partial K/\partial\beta\\. With the observed information, that
derivative is \\-\ell'''\\, which every family carries in closed form.
With the expected one it is
\\-\partial\\\mathbb{E}\[\ell''\]/\partial\eta\\, which is not
\\-\mathbb{E}\[\ell'''\]\\: differentiating an expectation moves the
measure as well as the integrand, so that
\\\partial\mathbb{E}\[\ell\_{ab}\]/\partial\eta_c =
\mathbb{E}\[\ell\_{abc}\] + \mathbb{E}\[\ell\_{ab}\ell\_{c}\]\\. The
mixed moment \\\mathbb{E}\[\ell\_{ab}\ell\_{c}\]\\ is not available from
a Bartlett identity on its own, which relates only symmetrized sums of
such moments.

distributions7 supplies it as
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html).
The two arrays are keyed differently, because the observed one is
symmetric in all three indices and the expected one in its first two
only, so the key builder travels with the array and a consumer does not
assume one key.

## See also

[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
