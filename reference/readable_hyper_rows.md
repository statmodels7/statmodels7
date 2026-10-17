# The Quantities That a Penalty's Hyperparameters Describe

Replaces the coordinate rows of a penalty whose hyperparameters are a
chart with the quantities that the penalty declares through
[`penalties7::penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.html):
the standard deviations and correlations of a correlated random effect,
in place of the logarithms of a Cholesky diagonal and the entries below
it.

## Usage

``` r
readable_hyper_rows(
  rd,
  th,
  Vh,
  p,
  key,
  level,
  role,
  src,
  cols,
  labels = character(0),
  at_edge = character(0)
)
```

## Arguments

- rd:

  The result of
  [`penalties7::penalty_readable()`](https://statmodels7.github.io/penalties7/reference/penalty_readable.html).

- th:

  The penalty's hyperparameters, as fitted.

- Vh:

  The hyperparameter variance matrix, or `NULL`.

- p:

  The distribution parameter whose equation contains the term.

- key:

  The penalty's key.

- level:

  The confidence level.

- role, src:

  The role and the source reported for the coordinate rows.

- cols:

  The column names of a summary block.

- labels:

  One label per coordinate of the matrix, or nothing. A correlation
  between the second and the third coordinate of a covariance is a
  number about two named effects, and the family that carries the chart
  does not identify which; where the labels are given, they replace the
  positions in the printed name.

- at_edge:

  The keys of the hyperparameter coordinates that
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
  found at a boundary. A quantity whose Jacobian entry on one of them is
  not structurally zero reports no standard error and no interval; every
  other quantity keeps its own.

## Value

A data frame of rows, in the shape of a summary block.

## Details

The standard error is the delta method, and it composes two Jacobians:
the penalty's, which is in the parameter scale of its hyperparameters,
and the link's, because the variance matrix is on the free scale on
which the outer criterion was maximized. Each interval is built on the
scale that the quantity declares, log for a standard deviation and
Fisher's z for a correlation, and mapped back, so a standard deviation
cannot be given a negative lower end and a correlation cannot be given
an interval that leaves \\(-1, 1)\\. The coefficients of a linear
predictor keep the symmetric interval of the estimate plus or minus
\\z\\ standard errors.

The summary does not print a test here, for the same reason that it does
not print one for the coordinate rows: a \\z\\ statistic tests the null
hypothesis that the quantity is zero, which for a standard deviation is
the edge of its range.
