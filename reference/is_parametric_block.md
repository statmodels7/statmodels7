# Whether a Term Is the Equation's Parametric Block

Returns `TRUE` for the block of unpenalized parametric columns, which is
the block that
[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md)
matches column by column.

## Usage

``` r
is_parametric_block(spec, param, nm)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- param:

  The distribution parameter.

- nm:

  The term's name in that parameter's equation.

## Value

A single logical.

## Details

The test reads the class of the term and not the name of the term. A
caller may give a term any name, and the deciding property is that the
block carries ordinary model-matrix columns whose names are variables
and levels, so that a coefficient of one model has the same meaning as
the coefficient of the same name in another. Every other kind of block
is a basis, a set of levels or a standardized design, where a shared
name does not imply a shared meaning.

## See also

[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
which uses it to determine how a block is matched.
