# Check the Argument That Selects a Variance Convention

Returns `x` if it is a
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
object and signals an error that lists the three constructors otherwise.

## Usage

``` r
as_vcov_method(x, arg = "vcov_type")
```

## Arguments

- x:

  The value of a `vcov_type` argument.

- arg:

  The argument's name, for the message.

## Value

`x`.
