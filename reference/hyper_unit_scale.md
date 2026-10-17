# The Units of the Free Values of a Matrix Prior

Computes one multiplier per row of `idx`: the factor that turns the
outer curvature in a hyperparameter's free value into the curvature in a
free value without units. It is 1 except for a free value of a matrix
prior that scales with the matrix.

## Usage

``` r
hyper_unit_scale(spec, design, hy, idx)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- hy:

  The fitted hyperparameters.

- idx:

  The hyperparameter index of
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

## Value

A numeric vector of length `nrow(idx)`.

## Details

A free value of a chart is homogeneous of degree \\k\\ in the scale of
the matrix when multiplying the matrix by \\\lambda^2\\ multiplies it by
\\\lambda^k\\. The off-diagonal free values of
[`parameters7::log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.html)
have \\k = 1\\ for a covariance, so they carry the units of the effects,
while a log-scale or a correlation coordinate has \\k = 0\\. The degree
is read from the chart, at a test point, with
[`parameters7::param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.html)
at \\\lambda = 2\\. The free value divided by \\c^k\\, with \\c^2\\ the
mean diagonal of the fitted matrix, has no units, and the curvature in
it is \\A\_{jj} (c^2)^k\\.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
which reads a boundary on the curvature multiplied by it.
