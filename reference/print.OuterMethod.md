# Print an Outer Method

Prints one line naming the criterion and the information it is built on,
as `<REML> observed information`, and for
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) the
number of folds and the selection rule.

## Usage

``` r
# S3 method for class 'OuterMethod'
print(x, ...)
```

## Arguments

- x:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- ...:

  Unused.

## Value

`x`, invisibly, as a print method should.

## See also

[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
and [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
for the objects printed.
