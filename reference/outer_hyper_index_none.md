# The Outer Index With No Hyperparameter

The empty index
[`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md)
returns where nothing is estimated, with the same columns and attributes
as a full one.

## Usage

``` r
outer_hyper_index_none()
```

## Value

A data frame of zero rows with columns `parameter`, `term` and `name`,
and the attributes `links` and `members`.

## Examples

``` r
nrow(statmodels7:::outer_hyper_index_none())
#> [1] 0
```
