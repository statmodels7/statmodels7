# Whether `par` Wrote Every Coordinate a Prior Drew

`TRUE` where every coordinate that a penalty's draw reached was then
written by `par`, so that the penalty's hyperparameters govern nothing
in the truth.

## Usage

``` r
targets_written(tg, ok, written)
```

## Arguments

- tg:

  The targets, from
  [`unit_draw_targets()`](https://statmodels7.github.io/statmodels7/reference/unit_draw_targets.md).

- ok:

  Which of them the draw reached.

- written:

  The result of
  [`par_written()`](https://statmodels7.github.io/statmodels7/reference/par_written.md).

## Value

A single logical.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md)
