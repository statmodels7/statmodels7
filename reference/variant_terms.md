# The Terms That Carry a Choice for a Search

Returns, for every term of the specification whose
[`modelterms7::term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.html)
is not `NULL`, the parameter, the call as written, the built term and
its pilot.

## Usage

``` r
variant_terms(spec)
```

## Arguments

- spec:

  A
  [StatmodSpec](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list, empty where no term needs a search.
