# The start_strategy Class Object

Returns the
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
class object from the namespace at the time of the call, so that a
membership test reads the class as it exists then.

## Usage

``` r
start_strategy_class()
```

## Value

The `start_strategy` S7 class object.

## Details

A class captured in a variable at build time makes such a test fail
under covr, which re-evaluates a package's code and so re-creates the
class.
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html)
against the captured copy is then `FALSE` for an object of the live
class. `distributions7` and `linkfunctions7` document the same behavior
of [`identical()`](https://rdrr.io/r/base/identical.html) on a class
that has been re-created.

## See also

[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
for the class itself.
