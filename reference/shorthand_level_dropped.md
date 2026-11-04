# The Level of a Shorthand gas() Beside an Intercept

Removes the development of the level from a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term whose subformulas came from the shorthand `by = ~f`, when the
equation's unpenalized columns span the constant.

## Usage

``` r
shorthand_level_dropped(term, others, n)
```

## Arguments

- term:

  An unbuilt term.

- others:

  The built additive terms of the same equation.

- n:

  The number of rows.

## Value

`term`, with the level's subformula removed where the rule above
applies.

## Details

The shorthand gives the same subformula to every parameter of the
filter. Beside an intercept, a development of the level \\\omega\\ would
duplicate it: its population value is confounded with the intercept and
would be held, and its departures per group would be a second set of
group levels, penalized through \\\omega\\ with a weight that falls as
the persistence approaches one. The level is therefore left undeveloped
and
[`statmod_held_levels()`](https://statmodels7.github.io/statmodels7/reference/statmod_held_levels.md)
holds it at zero; a level per cluster is a random intercept in the
equation. Without an intercept the development stays, and \\\omega\\
carries a population value and the departures. A development of
\\\omega\\ written as its own subformula is kept in either case.

## See also

[`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md),
the caller
