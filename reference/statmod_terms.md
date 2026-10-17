# Interpret and Build Each Parameter's Terms

Runs modelterms7's interpreter on every equation and builds the terms
that it names against the data.

## Usage

``` r
statmod_terms(equations, data, env, response = NULL, linpar = list())
```

## Arguments

- equations:

  A named list of one-sided formulas, one per distribution parameter.

- data:

  The data frame against which the terms are built.

- env:

  The environment that the original formula carried, which becomes the
  parent of the interpreting environment.

- response:

  The evaluated left-hand side, or `NULL` where there is none. Only
  [`seg_grid_start()`](https://statmodels7.github.io/statmodels7/reference/seg_grid_start.md)
  reads it.

- linpar:

  The options of the implicit parametric block, as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns them.

## Value

A list of two:

- `terms`:

  a named list with one entry per parameter, each a named list of built
  terms keyed by the term's call as written.

- `intercepts`:

  a named logical, one per parameter: whether that equation's parametric
  block carried an intercept.

## Details

The equations are interpreted with modelterms7's constructors in front
of the search path, so `s()` means this toolkit's whatever the caller
has attached. See
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
for the shim and the collision that it removes.

A factor covariate needs no handling of its own: the interpreter
collects the bare covariates of an equation into one
[`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html),
whose block comes from `model.matrix`, or from
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
when the block is sparse, and so carries the contrasts.

The starting positions of a break-point term are chosen on a grid,
through
[`seg_grid_start()`](https://statmodels7.github.io/statmodels7/reference/seg_grid_start.md),
when the caller did not name them, in place of the interior quantiles of
the covariate, which are the term's default.

## See also

[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
the caller,
[`terms_first()`](https://statmodels7.github.io/statmodels7/reference/terms_first.md)
for the shim.
