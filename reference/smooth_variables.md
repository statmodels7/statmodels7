# The Variables of a Smooth

Returns the names of the variables that a smooth term smooths over, in
the order of the call, without its `by` variable.

## Usage

``` r
smooth_variables(tm)
```

## Arguments

- tm:

  A
  [modelterms7::SmoothTerm](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.html).

## Value

A character vector with one name for each variable; a variable written
as an expression is named by its text.
