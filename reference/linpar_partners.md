# The Variables That Share a Linear Term with a Variable

Returns the variables of the terms of a linear formula that involve
`var`: `var` alone for a main effect, and the other variables of each
interaction in which it appears.

## Usage

``` r
linpar_partners(tm, var)
```

## Arguments

- tm:

  A
  [modelterms7::LinparTerm](https://statmodels7.github.io/modelterms7/reference/LinparTerm.html).

- var:

  A single string.

## Value

A character vector.
