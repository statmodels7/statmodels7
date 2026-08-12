# Which Shape of the Structural Contract a Term Implements

`"filter"` when the term reports a predictor, `"loglik"` when it reports
a likelihood, and `""` when it implements neither.

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

The question is asked of the methods the term registers, not of a list
of class names, so a structural term written outside the package is
routed without an edit here. The class a method was registered on is
`attr(m, "signature")[[1]]`, compared by name and package.
