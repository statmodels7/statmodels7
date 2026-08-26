# Print a Starting-Value Strategy

Prints a one-line description of a starting-value strategy: its label,
and for
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
whether the draw is centered, for
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
which optimizer it searches with.

## Arguments

- x:

  A
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
  or any object inheriting from it.

- ...:

  Unused.

## Value

`x`, invisibly, as a print method should.

## Examples

``` r
print(start_random())
#> <start> random around the intercept-only fit
```
