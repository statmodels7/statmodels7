# The Resolution That the Search Is Given, Out of the Readings Taken So Far

Returns the largest of the usable readings that
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
has returned during this search.

## Usage

``` r
resolution_summary(x)
```

## Arguments

- x:

  Every usable reading so far.

## Value

A single positive number, or `NA_real_` where there is none.

## Details

A reading is taken at every usable evaluation, and each one measures how
far the criterion moves when the mode is displaced by the error at which
that evaluation's inner fit stopped. The displacement varies across
evaluations by far more than the quantity that it represents does, so
the answer depends on which reading is kept.

The readings are already filtered to acceptable modes, because
[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md)
returns `NA` where the mode error exceeds
[`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md),
so every reading that reaches this function comes from a mode that the
layer accepts. The largest of them is therefore the one to keep. A
resolution has to bound how far the criterion can move between two such
modes, and not how little it moved at the evaluation whose inner score
happened to stop lowest.

The inner score at which a Newton step stops varies by orders of
magnitude across evaluations, and a minimum over such a sequence reads
the luckiest one, which understates the resolution. A resolution far
below the true one leaves the line search backtracking for improvements
that the criterion cannot resolve, and the run ends by exhausting
[`outer_backtracks()`](https://statmodels7.github.io/statmodels7/reference/outer_backtracks.md)
and not on a stopping rule.

## See also

[`criterion_resolution()`](https://statmodels7.github.io/statmodels7/reference/criterion_resolution.md),
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
