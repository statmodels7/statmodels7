# Check a Single Positive Number

Returns `x` invisibly when it is a single finite number above zero and
signals an error that names `arg` otherwise.

## Usage

``` r
check_positive_number(x, arg)
```

## Arguments

- x:

  The value to check.

- arg:

  The argument name used in the message.

## Value

`x`, invisibly.

## See also

[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md),
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`solve_qr()`](https://statmodels7.github.io/statmodels7/reference/solve_qr.md)
