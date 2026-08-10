# Interpret and Build Each Parameter's Terms

Runs modelterms7's interpreter on every equation and builds the terms it
names against the data.

## Usage

``` r
statmod_terms(equations, data, env)
```

## Arguments

- equations:

  A named list of one-sided formulas.

- data:

  A data frame.

- env:

  The environment the original formula carried.

## Value

A list with `terms` (a named list per parameter) and `intercepts` (a
named logical).

## Details

The equations are interpreted with modelterms7's constructors in front
of the search path, so that `s()` means ours whatever the user has
attached. A factor covariate needs no special handling: the interpreter
collects bare covariates into one `linpar()`, whose block comes from
`model.matrix` and therefore carries the contrasts.
