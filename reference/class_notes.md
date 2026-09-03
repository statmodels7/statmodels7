# What a Summary Says About the Covariance Classes

One note per class spanning more than one term, naming the label, the
grouping and the terms whose coefficients share the block.

## Usage

``` r
class_notes(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A character vector, possibly empty.

## Details

A class's hyperparameters are printed once, under its first member, so
without the note a reader sees a covariance of four coordinates under a
term carrying two columns and nothing saying where the other two came
from.

A class of one member gets no note: there is nothing shared to report,
and its block is the random effect it would have been without a label.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
which collects it.
