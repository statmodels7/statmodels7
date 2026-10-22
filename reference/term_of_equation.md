# A Term of an Equation, a Joint Filter Included

Returns the term kept under a key in the equation of a parameter, or the
joint filter of that key that spans the parameter where the filter is
kept in another equation.

## Usage

``` r
term_of_equation(spec, p, nm)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- p:

  The distribution parameter.

- nm:

  The key.

## Value

The term, or `NULL`.
