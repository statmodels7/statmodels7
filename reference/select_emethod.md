# The Method That a Criterion Reads

Returns the method of the expected information that the fits made for a
criterion use: the `method` of the first of `sparse_method` and `method`
whose information is an
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
object, and `emethod`, that of the inner method, when neither is.

## Usage

``` r
select_emethod(method, sparse_method, emethod)
```

## Arguments

- method, sparse_method:

  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  objects, or `NULL`.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

`NULL` or an object that inherits from
[distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).

## See also

[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
