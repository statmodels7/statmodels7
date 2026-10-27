# Unshrunk Estimates of the Effects of a Random Term

For each level, the least squares estimate of its coefficients from the
partial residuals of its observations, weighted by the prior weights.

## Usage

``` r
raw_effects(fit, param, tm, rb, information)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- param:

  The equation.

- tm:

  The random term.

- rb:

  The list from
  [`random_block()`](https://statmodels7.github.io/statmodels7/reference/random_block.md).

- information:

  An object from
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  or
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).

## Value

A matrix with the shape of `rb$b`; `NA` where a level cannot identify
its coefficients.
