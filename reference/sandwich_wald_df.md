# Degrees of Freedom of One Coefficient Under a CR2 Variance

Degrees of Freedom of One Coefficient Under a CR2 Variance

## Usage

``` r
sandwich_wald_df(fit, vcov_type, name)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vcov_type:

  The
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- name:

  The name of the coefficient, `"parameter:coefficient"`.

## Value

The Satterthwaite degrees of freedom when `vcov_type` is a frequentist
variance with
[`hc2_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md)
clustered and `reference = "t"`, and `Inf` otherwise.
