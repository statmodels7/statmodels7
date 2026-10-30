# The HC4 Correction of a Leverage Block With One Exponent per Equation

Computes \\F = \exp(-\tfrac12 D \log(I - H) D)\\ with \\D =
\mathrm{diag}(\delta^{1/2})\\, the correction of
[`hc4_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md).
It is \\\mathrm{diag}\\(1 - h\_{aa})^{-\delta_a/2}\\\\ when \\H\\ is
diagonal and \\(I - H)^{-\delta/2}\\ when the exponents are equal. Where
an eigenvalue of \\H\\ reaches one the logarithm does not exist, and the
block is corrected by
[`leverage_power()`](https://statmodels7.github.io/statmodels7/reference/leverage_power.md)
with the mean of the exponents.

## Usage

``` r
leverage_power_eq(H, delta)
```

## Arguments

- H:

  A symmetric leverage block.

- delta:

  The exponents of the equations.

## Value

A list with the matrix `F` and `capped`, the number of eigenvalues
treated as one.
