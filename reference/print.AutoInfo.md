# Print an Information Object Chosen Against the Family

Writes one line naming the object and the method of the expectation.

## Usage

``` r
# S3 method for class 'AutoInfo'
print(x, ...)
```

## Arguments

- x:

  An
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly. The function is called for the output that it writes.

## See also

[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)

## Examples

``` r
auto_info()
#> Information chosen against the family (closed form where the family has one, otherwise the outer product of the scores)
auto_info(method = distributions7::expect_mc(nsim = 500))
#> Information chosen against the family (expectation by Monte Carlo (nsim = 500))
```
