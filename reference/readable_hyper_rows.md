# The Quantities a Penalty's Hyperparameters Are About

Replaces the coordinate rows of a penalty whose hyperparameters are a
chart with the quantities it declares through
[`penalty_readable`](https://statmodels7.github.io/penalties7/reference/penalty_readable.html):
the standard deviations and correlations of a correlated random effect,
rather than the logarithms of a Cholesky diagonal and the entries below
it.

## Usage

``` r
readable_hyper_rows(rd, th, Vh, p, key, level, role, src, cols)
```

## Arguments

- rd:

  The result of
  [`penalty_readable`](https://statmodels7.github.io/penalties7/reference/penalty_readable.html).

- th:

  The penalty's hyperparameters, as fitted.

- Vh:

  The hyperparameter variance matrix, or `NULL`.

- p:

  The parameter the term sits in.

- key:

  The penalty's key.

- level:

  The confidence level.

- role, src:

  What the coordinate rows reported.

- cols:

  The column names of a summary block.

## Value

A data frame of rows, in the shape of a summary block.

## Details

The standard error is the delta method, and it composes two Jacobians:
the penalty's, which is in the parameter scale of its hyperparameters,
and the link's, the variance matrix being on the free scale the outer
criterion was maximized on. Each interval is built on the scale the
quantity declares – log for a standard deviation, Fisher's z for a
correlation – and mapped back, so a standard deviation cannot be given a
negative lower end and a correlation cannot be given an interval that
leaves \\(-1, 1)\\. That is the rule every other interval in the toolkit
follows.

No test is printed, for the reason the coordinate rows print none: the
null a \\z\\ of value over standard error reports on is that the
quantity is zero, which for a standard deviation is the edge of its
range.
