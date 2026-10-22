# The Joint Filters That Span an Equation

Returns the joint filters of a specification that filter the predictor
of a parameter, wherever the filter is kept, with the key under which
each is kept.

## Usage

``` r
joint_filters_of(spec, p)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- p:

  The distribution parameter.

## Value

A list of entries with `key`, `param` (the equation in which the filter
is kept) and `term`.
