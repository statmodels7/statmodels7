# The Coefficients That the Nested Expansion Predicts at a Step

Develops the unknown quantities of the nested problem as a polynomial in
the step of the held ones, and returns the point at its end: the
coefficients at their penalized mode, the hyperparameters, and the
specification that holds the coefficients that the criterion does not
integrate.

## Usage

``` r
profile_nested_predict(N, delta, order, se)
```

## Arguments

- N:

  The list that
  [`profile_nested_state()`](https://statmodels7.github.io/statmodels7/reference/profile_nested_state.md)
  returns.

- delta:

  The step of the held values, one number for each label of `which`.

- order:

  One more than the degree of the polynomial.

- se:

  The scale of each label of `which`.

## Value

A list with `par`, `hyper`, `spec` and `design`.

## Details

The coefficient of degree 1 is \\-V\_{uu}^{-1} V\_{uh}\delta\\. For a
higher degree \\j\\, the exact gradient of the criterion in the unknown
quantities is evaluated along the step, at the penalized mode of the
other coefficients, with the unknown quantities at the polynomial of
lower degree. Its \\j\\-th derivative at the origin is taken by one
central stencil of numericals7, and the coefficient of degree \\j\\ is
minus \\V\_{uu}^{-1}\\ times that derivative. Where nothing is unknown
the point is the penalized mode at the held values.
