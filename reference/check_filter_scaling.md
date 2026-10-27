# Check That a Family Supports a Scaled Score

Signals an error unless the family is univariate and its second
derivative of the expected information can be evaluated at the first
observation. The derivatives of a filter driven by a scaled score read
that second derivative, and the third and fourth derivatives of the
information are finite differences of it.

## Usage

``` r
check_filter_scaling(d, y, theta)
```

## Arguments

- d:

  A distribution object.

- y:

  The response.

- theta:

  The parameters, one value per observation or one in all.

## Value

`NULL`, invisibly.

## Details

The information that scales the score is part of the model, and it is
read in closed form at the predictor of every step of the recursion. A
family that carries only its density has no closed form, and neither of
the two substitutes defines the same model at a usable cost:

- [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  reads the information at the data as the square of the score,
  \\s_t^2\\, so the scaled score becomes \\s_t/\|s_t\|^{2d}\\, that is
  \\1/s_t\\ at \\d = 1\\ and the sign of \\s_t\\ at \\d = 1/2\\, which
  is another model;

- [`distributions7::expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.html)
  gives the information of the same model, but the recursion would
  evaluate one quadrature per step, and per derivative order, in
  sequence, since each step starts from the predictor of the step
  before.

A family supplies the closed form through methods for
[`distributions7::distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.html),
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html).
