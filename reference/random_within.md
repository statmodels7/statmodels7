# The Within-Group Rows of a Random-Effect Term

The row of the within-group design each observation carries, at the new
data where the specification has some and at the fitting rows otherwise.

## Usage

``` r
random_within(spec, param, key)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- param:

  The distribution parameter.

- key:

  The term's key.

## Value

A numeric matrix with one row per observation and one column per
within-group coordinate.
