# Resolve the Variance Convention Once

Returns `vt` unchanged unless it is a
[`bootstrap_vcov()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_vcov.md)
object, which is resolved by
[`boot_resolve()`](https://statmodels7.github.io/statmodels7/reference/boot_resolve.md).
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
call it on entry so that the replicates are computed once, however many
times they read the variance.

## Usage

``` r
resolve_vcov_method(fit, vt)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vt:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

## Value

A
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
object.
