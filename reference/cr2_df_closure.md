# Satterthwaite Degrees of Freedom of the CR2 Variance

Returns a function that computes the Satterthwaite degrees of freedom of
Pustejovsky and Tipton (2018) for the coefficients of the model, as
described at
[`hc2_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md).

## Usage

``` r
cr2_df_closure(parts, blk, Vb)
```

## Arguments

- parts:

  The list of
  [`leverage_parts()`](https://statmodels7.github.io/statmodels7/reference/leverage_parts.md).

- blk:

  The list of clusters built by
  [`sandwich_adjusted()`](https://statmodels7.github.io/statmodels7/reference/sandwich_adjusted.md),
  each with the rows `zr`, the design `Zg`, the roots `Whg` and the
  matrix \\(I - H\_{gg})^{-1/2}\\ as `F`.

- Vb:

  The inverse penalized information.

## Value

A function of `k`, positions among the coordinates of `Vb`, and
`contrast`, that returns a numeric vector of degrees of freedom. With
`contrast` (a matrix with one row for each coordinate of `Vb` and one
column for each quantity), the degrees of freedom are those of the
contrasts \\c'\beta\\ given by the columns, and `k` is ignored.
