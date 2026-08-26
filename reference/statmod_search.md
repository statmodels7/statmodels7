# How a Term Covers Its Own Hyperparameters

`"grid"` for every combination of the term's kinked hyperparameters,
`"cyclic"` for one at a time, or the default where the term named
neither.

## Usage

``` r
statmod_search(spec, row, default = "grid")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- row:

  One row of
  [`path_rows()`](https://statmodels7.github.io/statmodels7/reference/path_rows.md)'s
  index.

- default:

  What to use where the term named nothing.

## Value

`"grid"` or `"cyclic"`.

## Details

It belongs to the term, for the reason the whole enumeration does: a
criterion is asked about every hyperparameter of the model, and a smooth
one is read at the mode instead of being swept, so most of what it is
asked about could not answer. A penalty with a kink is fitted by a
scheme of its own, and how that scheme covers the term's own
hyperparameters is part of the scheme.

Being per term is also what keeps one term's choice off another's:
`y ~ lasso(X) + enet(R, search = "cyclic")` sweeps the elastic net one
coordinate at a time and leaves the lasso alone.

## See also

[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html),
[`statmod_path()`](https://statmodels7.github.io/statmodels7/reference/statmod_path.md)
