# The Reference Quantile of a Wald Interval

Returns the critical value for each coefficient of a table: the normal
quantile, or the \\t\\ quantile where the variance carries Satterthwaite
degrees of freedom.

## Usage

``` r
wald_reference(V, nm, se, level)
```

## Arguments

- V:

  A variance matrix, possibly with the attribute `"df_fn"` that
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
  copies from the sandwich.

- nm:

  The names of the coefficients, rows of `V`.

- se:

  Their standard errors; a coordinate without a finite one gets no
  degrees of freedom computed.

- level:

  The confidence level.

## Value

A list with `z` (a vector of critical values) and `df` (the degrees of
freedom, `Inf` where the reference is the normal).
