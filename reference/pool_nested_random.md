# A Term With Its Nested Random Effects Read in One Group's Columns

Returns a copy of a term carrying developed parameters in which every
row reads the within-group design of the listed random-effect sub-terms
in the columns of their first level, and nothing in the columns of the
other levels.

## Usage

``` r
pool_nested_random(tm, rows, coef)
```

## Arguments

- tm:

  A built term carrying developed parameters.

- rows:

  The rows of
  [`nested_random_terms()`](https://statmodels7.github.io/statmodels7/reference/nested_random_terms.md)
  for this term.

- coef:

  The term's coefficients to record as the ones it was last committed
  at, so that its block at new rows is the Jacobian there.

## Value

The modified term.

## Details

A random effect set aside contributes \\z_i^\top u\\ at every row, with
\\z_i\\ the row's within-group design and \\u\\ the value it is read at:
zero for `random = "zero"`, a node of the prior for `"marginal"`. Inside
a nonlinear term the effect enters the parameter, and the contribution
of the term is not linear in it, so the columns cannot simply be set to
zero as they are for a term written in an equation. With every row in
the first level, writing \\u\\ into that level's coefficients gives
every row the effect \\z_i^\top u\\, and the term is evaluated at those
coefficients by its own methods. A group the fit never saw is read the
same way, so it needs no effect of its own.

Two places are rewritten. At the fitting rows the term reads the stored
design of the developed parameter, `blueprint$Z`, whose columns for the
sub-term become the within-group rows (the sum of the group's columns, a
row being non-zero in its own group's alone). At new rows the term
reapplies each sub-term through
[`modelterms7::term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.html),
so the sub-term's grouping expression becomes the first level repeated.
