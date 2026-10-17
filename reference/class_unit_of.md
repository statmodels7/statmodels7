# The Covariance Class to Which One Term Belongs

Returns the penalized unit of the covariance class that covers a term's
coefficients, or `NULL` where the term has no label.

## Usage

``` r
class_unit_of(spec, design, param, nm)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- param:

  The distribution parameter in which the term sits.

- nm:

  The term's name in the formula.

## Value

One unit, as
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
returns it, or `NULL`.

## Details

A term with a label declares no penalty of its own, because the class
carries it, so the hyperparameters of such a term are read from the
class. A class's pieces record the pair of the equation and the term's
name, and the lookup uses that pair.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md)
for the classes,
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
for the reader that needs this.
