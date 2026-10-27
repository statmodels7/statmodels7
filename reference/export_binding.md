# The Active Binding of One Export of modelterms7

Returns a function of no arguments that gets the object called `name`
from the namespace of modelterms7. It is the function that
[`makeActiveBinding()`](https://rdrr.io/r/base/bindenv.html) takes in
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md).

## Usage

``` r
export_binding(name)
```

## Arguments

- name:

  A single string, the name of an export.

## Value

A function of no arguments.

## Details

The function is created here and not inside
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
so that its environment holds `name` and nothing else: a closure made in
the body of
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
would carry the environment that the caller passed.

## See also

[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
