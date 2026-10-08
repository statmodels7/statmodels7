# A Quadratic Form Row by Row, Where the Variance Has Missing Entries

\\x\_{ia}^\top V x\_{ib}\\ for every row \\i\\, with the entries of
\\V\\ that are not finite read as missing: a row is `NA` only where it
reaches one of them, that is where both its coefficient in \\x_a\\ and
its coefficient in \\x_b\\ are non-zero.

## Usage

``` r
row_quad(Xa, V, Xb)
```

## Arguments

- Xa, Xb:

  Matrices with a row per observation and a column per coefficient of
  the two blocks.

- V:

  The block of the variance matrix, rows for `Xa` and columns for `Xb`.

## Value

A numeric vector, one value a row.

## Details

A coefficient a kinked prior holds at its kink, such as a random effect
under a Laplace prior estimated at exactly zero, has no variance, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
reports `NA` there. A prediction whose design does not reach that
coefficient, the typical group's under `random = "zero"` for one, has a
standard error all the same; requiring every entry of the block to be
finite made it `NA` as well.
