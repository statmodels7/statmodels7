# S7 Class for the Meat of a HAC Variance

The class that
[`hac_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hac_sandwich.md)
instantiates. It carries the kernel, the bandwidth, and for a model
without a structural term the time and the group variables. Use the
constructor, which sets the properties.

## Usage

``` r
HacSandwich(
  type = character(0),
  kernel = character(0),
  bandwidth = NULL,
  time = NULL,
  cluster = NULL,
  data = NULL
)
```

## Arguments

- type:

  The string `"HAC"`.

- kernel:

  A single string, `"bartlett"`, `"parzen"` or `"truncated"`.

- bandwidth:

  `NULL` or a non-negative whole number.

- time, cluster:

  `NULL` or a one-sided formula with one variable.

- data:

  `NULL` or the data frame that holds those variables.

## Value

An S7 object with the properties `type`, `kernel`, `bandwidth`, `time`,
`cluster` and `data`.

## See also

[`hac_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hac_sandwich.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)

## Examples

``` r
hac_sandwich(bandwidth = 3)@kernel
#> [1] "bartlett"
```
