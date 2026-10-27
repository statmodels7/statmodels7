# The Information of the Leverage

Returns the information on which the hat matrix of a corrected sandwich
variance is computed: the expected information when the variance object
selects none, and the selected one otherwise.

## Usage

``` r
leverage_information(vcov_type, info)
```

## Arguments

- vcov_type:

  The
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- info:

  The information of the bread, as
  [`vcov_information()`](https://statmodels7.github.io/statmodels7/reference/vcov_information.md)
  resolves it.

## Value

An object from
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
or
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html).
