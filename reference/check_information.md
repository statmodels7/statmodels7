# Check That an Argument Is an Information Object

Returns `x` invisibly when it inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html)
and signals an error otherwise. With `auto = FALSE`, an
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
object is rejected as well, and with `null = TRUE`, `NULL` is accepted.

## Usage

``` r
check_information(
  x,
  arg = "information",
  fn = "reml",
  auto = TRUE,
  null = FALSE
)
```

## Arguments

- x:

  The value to check.

- arg:

  The argument name used in the message.

- fn:

  The function whose argument it is, used in the message.

- auto:

  Whether an
  [`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
  object is admissible.

- null:

  Whether `NULL` is admissible.

## Value

`x`, invisibly.

## See also

[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
