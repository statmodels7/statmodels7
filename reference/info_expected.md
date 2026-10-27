# Predicates on an Information Object

`info_expected()` returns `TRUE` for an
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
object and `info_auto()` for an
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
object. `info_label()` returns the word that names the matrix
(`"observed"`, `"expected"` or `"auto"`). `info_emethod()` returns the
method of the expected information that the object carries, `NULL` for
the closed form of the family or the outer product of the scores, and
`NULL` also for an
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
object.

## Usage

``` r
info_expected(x)

info_auto(x)

info_label(x)

info_emethod(x)
```

## Arguments

- x:

  An object inheriting from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html).

## Value

A logical, a string, `NULL` or an object inheriting from
[distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).

## See also

[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
