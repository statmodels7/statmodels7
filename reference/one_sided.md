# Build a One-Sided Formula From an Expression

Wraps a term expression as `~ expr` carrying the environment the
original formula had, so that a term's symbols resolve where the user
wrote them.

## Usage

``` r
one_sided(expr, env)
```

## Arguments

- expr:

  A language object.

- env:

  The environment to attach.

## Value

A one-sided formula.
