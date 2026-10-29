# The Coefficients That a Polynomial of the Trace Predicts at a Point

Develops the coefficients that are not held as a polynomial in the step
`delta` of the held ones, to degree `order - 1`, and returns its value
at the end of the step.

## Usage

``` r
profile_ray(S, delta, order)
```

## Arguments

- S:

  The list that
  [`profile_setup()`](https://statmodels7.github.io/statmodels7/reference/profile_setup.md)
  returns.

- delta:

  The step of the held coefficients, one number for each.

- order:

  One more than the degree of the polynomial.

## Value

The stacked coefficients at the end of the step.

## Details

The coefficient of degree 1 is \\-H\_{ff}^{-1} H\_{f\theta}\delta\\. For
each higher degree \\j\\, the gradient of the free coefficients is
evaluated along the step at the polynomial of degree \\j - 1\\, its
\\j\\-th derivative at the origin is taken by one central stencil of
numericals7 with a step that moves each held coefficient by a tenth of
its standard error, and the coefficient of degree \\j\\ is minus
\\H\_{ff}^{-1}\\ times that derivative.
