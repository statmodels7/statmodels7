# The Standard Deviation of a New Group's Parameter Under a Mixture

The Standard Deviation of a New Group's Parameter Under a Mixture

## Usage

``` r
group_sd_mixture(g, M, S, W, s)
```

## Arguments

- g:

  The link.

- M, S:

  The components' means and standard deviations of \\\eta\\, matrices
  with a row per observation and a column per component.

- W:

  The components' weights.

- s:

  The standard deviation of \\\eta^\*\\ itself, `NA` where it does not
  exist.

## Value

The standard deviation of \\h^{-1}(\eta^\*)\\, or `NA` where it is not
known to exist. See
[`group_interval()`](https://statmodels7.github.io/statmodels7/reference/group_interval.md).
