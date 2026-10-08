# Whether Fisher Scoring Is Contracting Slowly

`TRUE` where the last two ratios of consecutive scores both exceed
`rate`.

## Usage

``` r
iwls_slow(scores, rate)
```

## Arguments

- scores:

  The scores of the run so far, one per iteration.

- rate:

  The ratio above which a step counts as slow, or `NA` to never answer
  `TRUE`.

## Value

A single logical.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
[`iwls_switch_rate()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_rate.md).
