# Check That a Graph Has Its Data

Signals an error naming the graph when `data` is `NULL`, and when it
lacks one of `vars` or has a number of rows different from the fit.

## Usage

``` r
plot_needs_data(ctx, fit, type, vars = ctx$vars)
```

## Arguments

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- type:

  The
  [plot_type](https://statmodels7.github.io/statmodels7/reference/plot_type.md)
  object.

- vars:

  The variables the graph reads.

## Value

`ctx$data`, invisibly.
