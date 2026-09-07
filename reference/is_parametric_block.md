# Whether a Term Is the Equation's Parametric Block

`TRUE` for the block of unpenalized parametric columns, which is the one
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

The question is asked of the TERM rather than of its name. A caller may
name a term anything, and what decides is that the block carries
ordinary model-matrix columns whose names are variables and levels – so
a coefficient of one model means what the coefficient of the same name
means in another. Every other kind of block is a basis, a set of levels
or a standardized design, where a shared name is not a shared meaning.

## See also

[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md),
which uses it to decide how a block is matched.
