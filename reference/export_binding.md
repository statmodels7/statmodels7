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

The function is built from an expression, with the environment of the
base package, so that it holds nothing but the name. A closure written
in the body of
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
would carry the environment that the caller passed, and a closure
written here would carry its source reference: a package installed with
the source kept (`R_KEEP_PKG_SOURCE=yes`) serializes the source of the
whole package with each such reference.

## See also

[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
