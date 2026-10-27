# The Environment in Which the Time of a Term Is Evaluated

Returns the data of the graph as an environment's data when it was
given, and an empty list otherwise, so that the `time` expression of a
term is read from the data where it can be.

## Usage

``` r
plot_time_env(ctx, fit)
```

## Arguments

- ctx:

  The context described at
  [`plot_data()`](https://statmodels7.github.io/statmodels7/reference/plot_data.md).

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A list or a data frame.
