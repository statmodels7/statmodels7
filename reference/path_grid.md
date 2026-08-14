# The Values a Path Visits Over a Bounded Hyperparameter

An equally spaced grid strictly inside the hyperparameter's own
interval.

## Usage

``` r
path_grid(pen, name, n_values = 25L)
```

## Arguments

- pen:

  A penalties7 penalty.

- name:

  Which hyperparameter the path varies.

- n_values:

  How many points.

## Value

A numeric vector of values for `name`.

## Details

[`path_values`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
walks the SIZE OF THE KINK, from the value that empties the block down,
and a bounded hyperparameter cannot reach that end: the elastic net's
kink is \\\lambda\alpha\\, so at a given \\\lambda\\ no admissible
\\\alpha\\ empties the block and every point of such a path is dropped.
The interval itself is what a bounded shape is swept over instead –
\\\alpha\\ between the ridge and the lasso, SCAD's and MCP's shape
between their limits – and the sweeps being cyclic, one coordinate at a
time, the kink still moves through \\\lambda\\.

The endpoints are excluded because the bounds are open: the elastic net
at \\\alpha = 0\\ has no kink at all, and the path would be scoring a
penalty of another kind.

## See also

[`path_values`](https://statmodels7.github.io/statmodels7/reference/path_values.md),
[`path_bounded`](https://statmodels7.github.io/statmodels7/reference/path_bounded.md)
