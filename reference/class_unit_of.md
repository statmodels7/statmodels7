# The Covariance Class One Term Belongs To

The penalized unit of the covariance class that covers a term's
coefficients, or `NULL` where the term is not labelled.

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

  The distribution parameter the term sits in.

- nm:

  The term's name in the formula.

## Value

One unit, as
[`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md)
returns it, or `NULL`.

## Details

A labelled term declares no penalty of its own, the class carrying it,
so a reader asking a term for its hyperparameters has to ask the class
instead. The lookup is by the pair of the equation and the term's name,
which is what a class's pieces record.

## See also

[`statmod_classes()`](https://statmodels7.github.io/statmodels7/reference/statmod_classes.md)
for the classes,
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
for the reader that needs this.
