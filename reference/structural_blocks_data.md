# The Same Pieces as Data for the Compiled Recursion

Lays out the quantities that the callback of
[`.structural_blocks()`](https://statmodels7.github.io/statmodels7/reference/dot-structural_blocks.md)
reads as matrices, so that the compiled routes of `modelterms7` can read
them without calling back into R. The mixed second derivatives have one
column per distribution parameter (the filter's own column is zero and
the loop skips it). The third derivatives have one column per parameter
pair, and the pair \\(r, r_2)\\ is at column \\(r-1)\\n_p + r_2\\. The
static jacobian rows are densified, and the filter's parameter index is
included.

## Usage

``` r
structural_blocks_data(params, ap, Vs, H, D3, n, D4 = NULL, D5 = NULL)
```

## Arguments

- params:

  The distribution's parameter names.

- ap:

  Which of them carries the filter.

- Vs:

  The static rows.

- H, D3:

  The family's derivatives at the fitted predictors.

- n:

  The number of observations.

- D4, D5:

  The fourth and fifth derivatives, which the compiled third and fourth
  orders read: one column per index triple \\(r, r_2, r_3)\\ after the
  filter's own, at \\((r-1) n_p + r_2 - 1) n_p + r_3\\, and likewise for
  the quadruples of the fifth. `NULL` leaves them out.

## Value

A list with `H`, `D3`, `Vs` and `ap`, and `D4` and `D5` where they were
given.
