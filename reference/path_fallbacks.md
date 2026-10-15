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

The constructors
[`modelterms7::lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.html),
[`modelterms7::enet()`](https://statmodels7.github.io/modelterms7/reference/enet.html),
[`modelterms7::scad()`](https://statmodels7.github.io/modelterms7/reference/scad.html)
and
[`modelterms7::mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.html)
carry these settings on their own signatures, for example
`lasso(x, n_lambda = 25, min_ratio = 1e-4)` and
`enet(x, n_lambda = 25, n_alpha = 5)`, so that nothing here is used for
them. The settings here are used by a term that declares a kinked
penalty without an argument for the grid, such as
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
under a Laplace prior, whose hyperparameters are those of the effects'
distribution.

`kink` is the number of points on the axis swept by the size of the
kink, which runs geometrically over a factor of `1/min_ratio` (four
decades at the default). `other` is the number of points on an axis that
is swept over an interval instead, such as \\\alpha\\ between the ridge
and the lasso or the shape of SCAD or MCP over its usual range; with a
product of axes, every point added there multiplies the number of fits.

## See also

[`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md),
[`statmod_grid_size()`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md),
[`modelterms7::term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.html)
