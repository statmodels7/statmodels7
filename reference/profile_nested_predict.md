# The Coefficients That the Nested Expansion Predicts at a Step

Develops the coordinates that the criterion estimates as a polynomial in
the step of the held coefficients, and returns the coefficients of the
model at its end: those coordinates at the polynomial, the held ones at
the step, and the others at their penalized mode.

## Usage

``` r
profile_nested_predict(ctx, N, delta, order, se)
```

## Arguments

- ctx:

  The context of the call.

- N:

  The list that
  [`profile_nested_state()`](https://statmodels7.github.io/statmodels7/reference/profile_nested_state.md)
  returns.

- delta:

  The step of the held coefficients.

- order:

  One more than the degree of the polynomial.

- se:

  The standard errors of the held coefficients.

## Value

The stacked coefficients.

## Details

The coefficient of degree 1 is \\-V\_{\gamma\gamma}^{-1}
V\_{\gamma\theta} \delta\\. For a higher degree \\j\\, the exact
gradient of the criterion in the estimated coordinates is evaluated
along the step, at the penalized mode of the other coefficients, with
the coordinates at the polynomial of lower degree. Its \\j\\-th
derivative at the origin is taken by one central stencil of numericals7,
and the coefficient of degree \\j\\ is minus \\V\_{\gamma\gamma}^{-1}\\
times that derivative.
