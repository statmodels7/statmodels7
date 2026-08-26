# What a Path Does Where the Term Says Nothing

The number of values and the depth used for a kinked hyperparameter
whose term named neither.

## Usage

``` r
path_fallbacks()
```

## Value

A named list of `kink`, `other` and `min_ratio`.

## Details

The five penalized constructors carry these on their own signatures,
where a reader can see them, `lasso(x, n_lambda = 25, min_ratio = 1e-4)`
and `enet(x, n_lambda = 25, n_alpha = 5)`, so nothing here is reached
for them. What reaches it is a term that declares a kinked penalty
without offering an argument for the grid:
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
under a Laplace prior is the case, its hyperparameters being whatever
the effects' distribution happens to carry.

`kink` is the length of the path over the size of the kink, which runs
geometrically over `1/min_ratio`, which is four decades, and wants that
many points to be smooth in. `other` serves an axis that spans one
bounded interval instead, \\\alpha\\ between the ridge and the lasso or
a shape over its useful range, and needs fewer; with a product every
extra point there multiplies the fits.

## See also

[`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md),
[`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md),
[`modelterms7::term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.html)
