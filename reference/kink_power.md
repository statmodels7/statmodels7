# How the Size of the Kink Scales With a Hyperparameter

The exponent \\k\\ in \\s(v) = c\\v^k\\, read from the size of the kink
at the current value and at twice it, or `NA` where there is no such
exponent to read.

## Usage

``` r
kink_power(pen, theta, name)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one.

## Value

A list of the exponent `k`, the value `v0` it was read at and the size
`s0` there; `NULL` where the size does not move.

## Details

A separable penalty is minus a log density and its hyperparameter enters
as a scale, so the width of the subdifferential at the kink is a power
of it. Measured over four decades, the exponent is exactly one for the
lasso, for SCAD and MCP in \\\lambda\\ and for the elastic net in both
\\\lambda\\ and \\\alpha\\, and exactly zero for the shapes of SCAD and
MCP, which do not move the kink at all; the largest spread across a
decade sweep is \\5.6\times 10^{-16}\\. A Laplace prior written by its
scale has exponent \\-1\\, its kink narrowing as the hyperparameter
grows.

The exponent is measured, never assumed, and whoever inverts it checks
the answer, so a penalty whose kink is not a power of its hyperparameter
costs two evaluations and then takes the search.

## See also

[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
