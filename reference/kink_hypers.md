# Which Hyperparameters Set the Size of the Kink

The names of the hyperparameters whose value changes
[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md).

## Usage

``` r
kink_hypers(pen, theta, unbounded = TRUE)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- unbounded:

  Whether to keep only those with an infinite upper bound. `TRUE` by
  default.

## Value

A character vector, possibly empty.

## Details

Each hyperparameter is moved to another admissible value and the size of
the kink is read again, so that no list of families is consulted. The
shape parameters of SCAD and MCP leave the subdifferential at zero
unchanged and govern how fast the penalty flattens further out, while
\\\lambda\\ and the elastic net's \\\alpha\\ both scale it.

`unbounded` restricts the answer to the hyperparameters with no upper
bound.
[`path_by_kink()`](https://statmodels7.github.io/statmodels7/reference/path_by_kink.md)
calls this function with `unbounded = FALSE` and excludes the bounded
hyperparameters itself; those are swept over their interval by
[`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md).

## See also

[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)
