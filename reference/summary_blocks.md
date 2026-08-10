# The Blocks of One Distribution Parameter

Groups a parameter's terms into the readings a summary prints: the
parametric terms together, and one block per penalized term.

## Usage

``` r
summary_blocks(fit, spec, design, p, ci)
```

## Arguments

- fit:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec:

  The specification.

- design:

  The design.

- p:

  The distribution parameter.

- ci:

  The flat interval table, as
  [`confint.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
  returns it with the statistic and the p-value added.

## Value

A list of block records, each with `kind`, `label`, `n_coef`, `edf`,
`n_zero` and `table`.
