# Gauss-Hermite Nodes and Weights

Returns the nodes and weights of the rule for \\\int f(x) e^{-x^2} dx\\,
computed by the eigendecomposition of the Jacobi matrix of the Hermite
polynomials (Golub and Welsch 1969).

## Usage

``` r
gauss_hermite(k)
```

## Arguments

- k:

  The number of nodes.

## Value

A list with `x` and `w`.
