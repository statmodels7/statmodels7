# How the Expected Information Moves With the Coefficients, Twice

\\\partial^2\\\mathbb{E}\[\ell''\]/\partial\eta^2\\, the array the
expected route's outer Hessian traces twice-contracted against the
leverage diagonal, together with its key builder.

## Usage

``` r
ctx_kmove2(ctx, spec, design, coef, hyper)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

## Value

A list with `deriv` and `key`.

## Details

It is read only on the expected route, where
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
has already asked
[`expected_deriv2_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv2_ok.md).
The components are symmetric in the information's pair and in the pair
differentiated in, separately, so they are read through
[`distributions7::d2expected_key()`](https://statmodels7.github.io/distributions7/reference/d2expected_key.html)
and never through the observed route's sorted quadruple.

## See also

[`ctx_kmove()`](https://statmodels7.github.io/statmodels7/reference/ctx_kmove.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
