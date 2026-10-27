# The Information Bound to an Inner Method

Returns the information object that
[`inner_bind()`](https://statmodels7.github.io/statmodels7/reference/inner_bind.md)
stored in `method`: the property `hessian` of an
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
object, and the attribute `"inner_information"` of an optimizers7
optimizer. An optimizer that was never bound has the observed
information.

## Usage

``` r
inner_info(method)
```

## Arguments

- method:

  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

## Value

An object that inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## See also

[`inner_bind()`](https://statmodels7.github.io/statmodels7/reference/inner_bind.md)
