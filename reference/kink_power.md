# How the Size of the Kink Scales With a Hyperparameter

Returns the exponent \\k\\ in \\s(v) = c\\v^k\\, read from the size of
the kink at the current value and at twice it (or at half it, where
twice it is outside the hyperparameter's interval), or `NULL` where
there is no such exponent to read.

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

A list of the exponent `k`, the value `v0` at which it was read and the
size `s0` there; `NULL` where the size does not move or where either
size or the value is not positive and finite.

## Details

A separable penalty is minus a log density and its hyperparameter enters
as a scale, so the width of the subdifferential at the kink is a power
of it. The exponent is one for the lasso, for SCAD and MCP in
\\\lambda\\ and for the elastic net in both \\\lambda\\ and \\\alpha\\.
For the shapes of SCAD and MCP, which do not move the kink, there is no
exponent and the result is `NULL`. A Laplace prior written by its scale
has exponent \\-1\\, and its kink narrows as the hyperparameter grows.

The exponent is measured and not assumed, and
[`kink_by_power()`](https://statmodels7.github.io/statmodels7/reference/kink_by_power.md)
checks the inversion, so that a penalty whose kink is not a power of its
hyperparameter costs two evaluations of the size and then takes the
search of
[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md).

## See also

[`kink_solve()`](https://statmodels7.github.io/statmodels7/reference/kink_solve.md),
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
