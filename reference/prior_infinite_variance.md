# Whether a Prior Has No Finite Variance

Whether a Prior Has No Finite Variance

## Usage

``` r
prior_infinite_variance(ob)
```

## Arguments

- ob:

  An entry of `other` from
  [`prior_parts()`](https://statmodels7.github.io/statmodels7/reference/prior_parts.md).

## Value

`TRUE` where the prior is univariate and the variance its family reports
is infinite or not a number (a Cauchy); `FALSE` otherwise, including
where no variance can be read.
