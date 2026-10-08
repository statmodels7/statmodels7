# The Contraction Rate Above Which Fisher Scoring Yields

The ratio of consecutive scores above which
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
continues on the observed information before
[`iwls_switch_after()`](https://statmodels7.github.io/statmodels7/reference/iwls_switch_after.md)
steps.

## Usage

``` r
iwls_switch_rate()
```

## Value

A single number, or `NA` for the count alone.

## See also

[`iwls_slow()`](https://statmodels7.github.io/statmodels7/reference/iwls_slow.md).
