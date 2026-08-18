# The Same Pieces as Data for the Compiled Recursion

The quantities the callback of
[`.structural_blocks()`](https://statmodels7.github.io/statmodels7/reference/dot-structural_blocks.md)
reads, laid out as matrices so `modelterms7`'s compiled second-order
route can read them without calling back into R: the mixed second
derivatives one column per distribution parameter (the filter's own
column zero, the loop skips it), the third derivatives one column per
parameter pair with pair \\(r, r_2)\\ at column \\(r-1)\\n_p + r_2\\,
the static jacobian rows densified, and the filter's parameter index.

## Usage

``` r
structural_blocks_data(params, ap, Vs, H, D3, n)
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

## Value

A list with `H`, `D3`, `Vs` and `ap`.
