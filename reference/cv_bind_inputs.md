# Carry a Term's Matrix Input Onto a Subset of the Rows

Adds, as a column of the subset, each matrix that a term was given, so
that rebuilding the model on those rows finds it there.

## Usage

``` r
cv_bind_inputs(spec, sub, i, n)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  whose terms are built.

- sub:

  The subset of the data, already taken.

- i:

  The rows used to take it.

- n:

  The number of rows on which the fit was built.

## Value

`sub`, which has a column per matrix input that the terms carry.

## Details

`interpret_formula()` evaluates a term's call as
`eval(call, data, env)`, so a name is looked up in `data` first and in
the formula's environment after. `data.frame(X = X, y = y)` splits a
matrix into `X.x1 ... X.xp`, leaving no column `X`, so `lasso(X)`
reaches past the data to the matrix in the calling environment. The fit
itself is correct (the matrix is captured once, and the coefficients are
identical to those of the other spelling), but the fold cannot evaluate
the term on its rows: the name still resolves to all the rows.

The matrix is already on the built term, and `term_build()` checked at
the full fit that it has one row per observation, so the rows of a fold
are the same rows by position. Binding the subset here gives the fold
the spelling that the documentation of the matrix input recommends, with
the matrix as a column of the data.

A matrix carries no knots, contrasts or levels, so subsetting it and
evaluating it again give the same block. A formula input is left
untouched.

It applies where `input_expr` is a plain symbol, which is the case for
which the name can be bound. A call such as `lasso(scale(X))` keeps only
its own value on the term, never the `X` that its re-evaluation would
need, so it is left to the error that names it.

## See also

[`cv_curve()`](https://statmodels7.github.io/statmodels7/reference/cv_curve.md)
