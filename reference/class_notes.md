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

A class's hyperparameters are printed once, at the head of the summary,
where each coordinate is named for the equation and the column it
belongs to. The note says the same thing in one sentence, for a reader
who has the summary object rather than the printed page.

A class of one member gets no note: there is nothing shared to report,
and its block is the random effect it would have been without a label.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md),
which collects it;
[`summary_class_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_class_blocks.md),
which prints the numbers.
