# The Numeric by Variable of a Smooth

Returns the name of the numeric `by` variable of a built
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::scs()`](https://statmodels7.github.io/modelterms7/reference/scs.html)
or
[`modelterms7::scte()`](https://statmodels7.github.io/modelterms7/reference/scte.html)
term, whose contribution is \\z\\f(x)\\, and an empty vector for any
other term. The graphs of such a term draw the coefficient \\f(x)\\, the
contribution at \\z = 1\\.

## Usage

``` r
numeric_by_variable(tm)
```

## Arguments

- tm:

  A model term.

## Value

A character vector of length 0 or 1.
