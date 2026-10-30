# The Common Layers of a Profile Panel

Draws the profile of one label against its value, with the estimate, the
limit of the interval at level 0.95 where the profile gives a likelihood
ratio statistic, and the points whose refit did not converge. Every
element is mapped to an aesthetic with a descriptive label, so that the
legend names it;
[`profile_finish()`](https://statmodels7.github.io/statmodels7/reference/profile_finish.md)
adds the legend.

## Usage

``` r
profile_base(
  q,
  x,
  column = x@profile_column[[q@label]],
  cutoff = x@is_likelihood_ratio[[q@label]],
  steps = FALSE,
  line = TRUE
)
```

## Arguments

- q:

  A
  [profiled_quantity](https://statmodels7.github.io/statmodels7/reference/profiled_quantity.md)
  object.

- x:

  The
  [StatmodProfile](https://statmodels7.github.io/statmodels7/reference/StatmodProfile.md)
  object.

- column:

  The column of `points` that is drawn.

- cutoff:

  Whether to draw the line of the interval at level 0.95.

- steps:

  Whether the profile is constant between two consecutive values of the
  covariate, and is drawn as one segment over each interval.

- line:

  Whether the points of the profile are joined by a line.

## Value

A list with `plot`, a ggplot2 object, `data`, the data frame drawn, and
`entries`, the legend entries of the colour scale, a list of lists with
`label`, `colour`, `shape` and `linetype`.
