# Draws from a Multivariate Normal Distribution

Draws `n` vectors from \\N(0, V)\\ through the eigendecomposition of
\\V\\, which also covers a singular \\V\\.

## Usage

``` r
mvn_draws(V, n)
```

## Arguments

- V:

  A symmetric positive semi-definite matrix.

- n:

  The number of draws.

## Value

A matrix with `n` rows and `ncol(V)` columns.
