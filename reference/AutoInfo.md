# S7 Class for the Information Chosen Against the Family

`AutoInfo` is the class of the object that
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
returns: an information whose matrix, observed or expected, is settled
once the family is known. It inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html),
the abstract parent of
[distributions7::ObservedInfo](https://statmodels7.github.io/distributions7/reference/Information.html)
and
[distributions7::ExpectedInfo](https://statmodels7.github.io/distributions7/reference/Information.html).
Use
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md),
which validates its arguments; the class is exported so that
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
can test an object.
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
compares the three objects.

## Usage

``` r
AutoInfo(method = NULL, repair_floor = integer(0))
```

## Arguments

- method:

  `NULL` or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html),
  the method that obtains the expected information wherever the settled
  choice or the fallback uses it. `NULL` selects the closed form of the
  family where it has one and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise.

- repair_floor:

  A single positive number, the relative floor of the repair of an
  indefinite observed information.

## Value

An S7 object of class `AutoInfo`, which inherits from
[distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html),
with the properties `method` and `repair_floor`.

## See also

[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
for the constructor and
[information_objects](https://statmodels7.github.io/statmodels7/reference/information_objects.md)
for the comparison with
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
and
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html).

## Examples

``` r
a <- auto_info(method = distributions7::expect_mc(nsim = 500))
S7::S7_inherits(a, distributions7::Information)
#> [1] TRUE
S7::S7_inherits(a, AutoInfo)
#> [1] TRUE
a@method
#> Expectation by Monte Carlo (nsim = 500)
```
