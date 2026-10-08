# The Coefficients That Move Along a Ridge

The rows of the coefficient table that lose their standard error to the
ridges
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
names.

## Usage

``` r
ridge_coef_rows(ridge, V, rows)
```

## Arguments

- ridge:

  The certificate's `ridge`, each element carrying the attribute
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
coefficient with \\\|b_j\| \ge 10^{-2}\\ is taken too, the direction
being scaled so that its largest entry is one. Where the variance matrix
has no finite entries for the ridge's coordinates (the mode reading
already held one of them), the named coefficients alone are taken.

## See also

[`certificate_ridges()`](https://statmodels7.github.io/statmodels7/reference/certificate_ridges.md),
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
