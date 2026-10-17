# The Properties That Every Criterion Carries

Returns the defaults for `nfolds`, `rule` and `folds`, which only
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
reads, so that one class can serve every criterion and each constructor
names only what it needs.

## Usage

``` r
outer_path_defaults()
```

## Value

A named list with `nfolds`, `rule` and `folds`, to be spliced into an
[`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
call.

## Details

The settings of a path are not among them. How many values a path
visits, how far down it reaches, and whether a term's own
hyperparameters are combined or swept one at a time all belong to the
term.
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
are also used for a model's smooth hyperparameters, which are read at
the mode and not swept, so a path setting on the criterion would not
apply to them.
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) is
accepted only as a `sparse_criterion`.

## See also

[`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
for the settings to which a term falls back when it declares a kinked
penalty without offering a grid argument,
[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html)
for where a path's shape is set.
