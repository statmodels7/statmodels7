# S7 Classes for the Shipped Test Statistics

These are the classes that
[`wald_test()`](https://statmodels7.github.io/statmodels7/reference/wald_test.md),
[`lr_test()`](https://statmodels7.github.io/statmodels7/reference/lr_test.md),
[`score_test()`](https://statmodels7.github.io/statmodels7/reference/score_test.md)
and
[`gradient_test()`](https://statmodels7.github.io/statmodels7/reference/gradient_test.md)
instantiate. Each inherits from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
and carries its `name` and `label`. `RestrictedTest` is the abstract
parent of the three classes whose statistic is evaluated at a restricted
refit (`LrTest`, `ScoreTest` and `GradientTest`), and
[`test_restricted()`](https://statmodels7.github.io/statmodels7/reference/test_restricted.md)
returns `TRUE` for every object that inherits from it. Use the
constructors, which set the properties.

## Usage

``` r
WaldTest(name = character(0), label = character(0))

RestrictedTest(name = character(0), label = character(0))

LrTest(name = character(0), label = character(0))

ScoreTest(name = character(0), label = character(0))

GradientTest(name = character(0), label = character(0))
```

## Arguments

- name, label:

  Single strings: the short name of the statistic and the name under
  which a result prints it.

## Value

An S7 object inheriting from
[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md),
with the properties `name` and `label`.

## See also

[test_statistic](https://statmodels7.github.io/statmodels7/reference/test_statistic.md)
for the comparison of the four statistics.

## Examples

``` r
S7::S7_inherits(lr_test(), test_statistic)
#> [1] TRUE
class(lr_test())
#> [1] "statmodels7::LrTest"         "statmodels7::RestrictedTest"
#> [3] "statmodels7::test_statistic" "S7_object"                  
```
