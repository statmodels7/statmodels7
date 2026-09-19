# A Design Times Coefficients, Over the Nonzero Ones

Computes \\X\beta\\ reading only the columns whose coefficient is not
zero, which after a lasso is most of a block.

## Usage

``` r
x_times_b(X, b)
```

## Arguments

- X:

  A design block, dense or a Matrix class.

- b:

  Its coefficients, `ncol(X)` numbers.

## Value

A numeric vector of `nrow(X)` entries.

## Details

A term \\0 \cdot x\_{ij}\\ adds exactly zero to the running sum, so
leaving those columns out gives the same number, the remaining products
being summed in the same column order. The design is subset only where
fewer than all the coefficients are nonzero; with none, the result is a
vector of zeros. A non-finite coefficient keeps the full product. The
design itself is taken to be finite, as a model matrix is: an infinite
entry in a column whose coefficient is zero would make the full product
`NaN` and this one not.

## See also

[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md),
[`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md)
