# The Coefficients That Move Along a Ridge

Returns the rows of the coefficient table whose standard errors are
dropped because of the ridges that
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
names.

## Usage

``` r
ridge_coef_rows(ridge, V, rows)
```

## Arguments

- ridge:

  The certificate's `ridge`, in which each element carries the attribute
  `where`, the stacked positions of its coefficients (`NA` for a
  hyperparameter).

- V:

  The variance matrix over the coefficients, or `NULL`.

- rows:

  The row names of the coefficient table, in stacked order.

## Value

The row indices, sorted.

## Details

A ridge is written over the coordinates of the outer criterion; the
coefficients among them are taken as they are. The coefficients held in
the inner fit move with the ridge through the mode, and the variance
matrix carries that: the change of coefficient \\j\\ per unit along the
combination \\x'\gamma\\ is \\b_j = (Vx)\_j / (x'Vx)\\, and a
coefficient with \\\|b_j\| \ge 10^{-2}\\ is taken too; the direction is
scaled so that its largest entry is one. Where the variance matrix lacks
a row for one of the ridge's coefficients or has a non-finite entry
among them (the mode reading already held one of them), only the
coefficients that the ridge names are taken.

## See also

[`certificate_ridges()`](https://statmodels7.github.io/statmodels7/reference/certificate_ridges.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
