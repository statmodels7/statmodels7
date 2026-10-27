# The Variables of a Term

Returns the names of the variables that a term reads: those of the
formula for the linear part, and those of the call that names the term
otherwise; for a nonlinear term, only those of the call that the build
recorded as data, since the call also names the term's parameters.

## Usage

``` r
term_variables(tm, name)
```

## Arguments

- tm:

  A model term.

- name:

  The name of the term in its equation.

## Value

A character vector.
