# S7 Class for the Meat of a Sandwich Variance

The class that
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
and
[`hc1_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
instantiate. It carries the type of the finite-sample correction, the
optional cluster variable and the optional data frame from which that
variable is read. Use the constructors, which set the properties.

## Usage

``` r
SandwichMethod(type = character(0), cluster = NULL, data = NULL)
```

## Arguments

- type:

  A single string, `"HC0"` or `"HC1"`.

- cluster:

  `NULL` or a one-sided formula with one variable, such as `~ id`.

- data:

  `NULL` or the data frame that holds the cluster variable.

## Value

An S7 object with the properties `type`, `cluster` and `data`.

## See also

[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)

## Examples

``` r
class(hc1_sandwich(cluster = ~ id))
#> [1] "statmodels7::SandwichMethod" "S7_object"                  
hc1_sandwich(cluster = ~ id)@type
#> [1] "HC1"
```
