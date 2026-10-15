# Does a Penalty Supply What a Marginal Criterion Needs?

Evaluates the derivatives of a penalty at a probe value of its
hyperparameters and reports whether they could be computed.

## Usage

``` r
penalty_answers(pen, order = 1L)
```

## Arguments

- pen:

  A penalties7 penalty.

- order:

  `1` or `2`.

## Value

A single logical.

## Details

The gradient requires
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html).
A penalty that is not quadratic in the coefficients
([`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html))
must also supply
[`penalties7::penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.html),
because its Hessian moves with the mode. The order-2 route additionally
requires
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html)
and
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)
and, for a penalty that is not quadratic,
[`penalties7::penalty_d2hessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.html)
and
[`penalties7::penalty_dhessian_beta_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian_beta.html).
A penalty without hyperparameters returns `FALSE`.
