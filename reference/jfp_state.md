# The State of the Jeffreys Penalty at One Point

Computes, at given coefficients, everything that the value, the gradient
and the Hessian of the Jeffreys penalty \\\sum_b a_b\log\|\tilde I_b\|\\
read: for each block its rows, its expected information, the Cholesky
factor of that information, the whitened rows and the leverages.

## Usage

``` r
jfp_state(spec, design, coef, blocks, emethod = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- blocks:

  A list of
  [`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
  objects.

- emethod:

  The `method` of
  [`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html),
  or `NULL` for the family's own closed form.

## Value

An environment with the fields `spec`, `design` (refreshed at `coef`),
`coef`, `ev`, `params`, `npar`, `offs`, `nb` (the number of
coefficients), `total` (the length of \\u\\), `joint` (`NULL`, or what
[`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md)
returns), `E`, `emethod` and `blocks`, each entry of the last carrying
`block`, `index` (positions in \\u\\), `rows`, `info`, `U`, `Z`, `G`,
`logdet` and `ok`.

## Details

For a block \\b\\ with rows \\\tilde V\_{q,t}\\, \$\$\tilde I_b = \sum_t
w_t\sum\_{q,r}\mathcal I\_{qr}(\eta_t)\\ \tilde V\_{q,t}^\top\tilde
V\_{r,t},\$\$ where \\\mathcal I\_{qr}\\ is the family's expected
information on the link scale, minus
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html).
With \\\tilde I_b = U^\top U\\, the whitened rows are \\Z_q = \tilde V_q
U^{-1}\\ and the leverages are \\h\_{qr,t} = Z\_{q,t}Z\_{r,t}^\top\\.

The point is the joint vector \\u\\ of
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md).
For a model without a structural term it is the stacked coefficients.
For a model carrying a filter it is the coefficients followed by the
filter's free parameters, and the derivative of every predictor in \\u\\
is read from the forward recursion
([`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md));
a block in an equation that the filter drives moves with the
coefficients and is rejected in this phase. A model carrying a mixture
over latent states is always rejected, because its expected information
does not exist.

The derivatives of the expected information are computed at the first
request and kept in the state, so that the value, the gradient and the
Hessian at one point call the family once per order.

Where a block's information is not positive definite its entry carries
`ok = FALSE` and a log-determinant of `-Inf`, and the gradient and the
Hessian signal an error.

## See also

[`jfp_value()`](https://statmodels7.github.io/statmodels7/reference/jfp_value.md),
[`jfp_gradient()`](https://statmodels7.github.io/statmodels7/reference/jfp_gradient.md),
[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
