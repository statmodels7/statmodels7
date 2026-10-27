# A Sum of Terms Drawn as One Effect

The class of the object on which
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md)
dispatches when a variable names more than one term of an equation, as
`x` does in `y ~ x + s(x)`. Each element of `parts` is a list with the
elements `name`, `term` and `cols`, as in the context of
[`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

## Usage

``` r
TermSum(parts = list())
```

## Arguments

- parts:

  A list of the parts.

## Value

A `TermSum` object.
