# Resolve the Term of a Graph

Finds the term, or the terms, that `term` names in the equation `param`,
and the columns of the linear part that belong to a variable.

## Usage

``` r
resolve_plot_term(fit, term, param, data)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- term:

  A single string, the name of a term or of a variable.

- param:

  The equation, or `NULL`.

- data:

  The data frame of the fit, or `NULL`; the columns of the linear part
  that belong to a variable are read from it.

## Value

A list with the elements `param`, `vars`, `title` and `parts`; each
element of `parts` is a list with `name`, `term` and `cols`.
