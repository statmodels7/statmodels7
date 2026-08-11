# Choose a Point of the Path

The index of the smallest criterion, or of the largest kink whose
criterion is within one standard error of it.

## Usage

``` r
path_pick(value, se = NULL, rule = "min")
```

## Arguments

- value:

  The criterion at each point.

- se:

  Its standard error, or `NULL`.

- rule:

  `"min"` or `"1se"`.

## Value

A single index, or `NA` where no point was usable.

## Details

The one-standard-error rule takes the sparsest fit that is not
measurably worse than the best one, which is Breiman's rule as glmnet
applies it. The path runs from the emptiest fit to the fullest, so the
largest kink among the admissible points is the first of them.

## References

Breiman, L., Friedman, J. H., Olshen, R. A. and Stone, C. J. (1984).
*Classification and Regression Trees*. Wadsworth.
