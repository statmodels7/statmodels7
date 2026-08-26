# The Properties Every Criterion Carries

The defaults for the three properties only
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

What a **path** does is not among them. How many values it visits, how
far down it reaches, and whether a term's own hyperparameters are
combined or swept one at a time all belong to the term. The same
criterion object is put to the model's smooth hyperparameters as well,
and those are read at the mode instead of being swept, so a path setting
on the criterion would be meaningless for most of what it is asked
about.

## See also

[`path_fallbacks()`](https://statmodels7.github.io/statmodels7/reference/path_fallbacks.md)
for the settings a term that declares a kinked penalty without offering
a grid argument falls back to,
[`modelterms7::term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.html)
for where a path's shape is set.
