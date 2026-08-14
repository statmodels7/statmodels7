# How Many Points an Axis Beside the Path Gets

The default number of values for a hyperparameter that is NOT swept by
the size of its kink, which is one fifth of `n_values` and at least two.

## Usage

``` r
path_n_other(method)
```

## Arguments

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

A single integer.

## Details

`n_values` is the length of the path over the size of the kink, which
runs geometrically over `1/min_ratio` – four decades at the defaults –
and wants that many points to be smooth in. An axis beside it spans one
bounded interval, \\\alpha\\ between the ridge and the lasso or a shape
over its useful range, and does not.

The ratio is a DEFAULT and not a derivation: with a product every extra
point on this axis multiplies the fits, so one fifth is what keeps
`search = "grid"` affordable at the shipped `n_values`, 25 by 5 rather
than 25 by 25. A term that wants otherwise says so –
`enet(x, n_alpha = 12)`, `scad(x, n_a = 8)` – and is obeyed.

## See also

[`path_grid`](https://statmodels7.github.io/statmodels7/reference/path_grid.md),
[`statmod_grid_size`](https://statmodels7.github.io/statmodels7/reference/statmod_grid_size.md)
