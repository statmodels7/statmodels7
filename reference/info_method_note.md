# The Method of an Information Object in Words

Returns the text that follows "Expected information" in a header: the
empty string for the closed form of the family (`method = NULL`), and
the method in parentheses otherwise.

## Usage

``` r
info_method_note(x)
```

## Arguments

- x:

  An
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
  object.

## Value

A single string.

## See also

[`info_title()`](https://statmodels7.github.io/statmodels7/reference/info_title.md)
