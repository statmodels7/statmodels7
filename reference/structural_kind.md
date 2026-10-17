# Which Shape a Structural Term Implements

Returns `"filter"` when the term registers its own
[`modelterms7::term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.html)
method, `"loglik"` when it registers its own
[`modelterms7::term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.html)
method, and `""` when it registers neither.

## Usage

``` r
structural_kind(term)
```

## Arguments

- term:

  A term.

## Value

A single string.

## Details

The shape is read from the methods that the term registers, not from a
list of class names, so a structural term written outside the package is
routed without an edit here. The class on which a method was registered
is `attr(m, "signature")[[1]]`, and a method registered on
`structural_term` itself does not count. The filter shape is tested
first.
