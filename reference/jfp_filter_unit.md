# A Jeffreys Block in the Subformula of a Filter's Parameter

Resolves a Jeffreys unit that a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term carries up from the subformula of one of its parameters: the
positions among the term's parameters from which the block's rows are
combined, and whether the constant direction of the subformula enters
them.

## Usage

``` r
jfp_filter_unit(u, spec)
```

## Arguments

- u:

  The unit, structural, with `within` the block's positions among the
  term's parameters.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A list with `kind` (`"filter"`), `const`, `basis` (positions among the
term's parameters) and `C`, as
[`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
takes them.

## Details

The rows are the derivatives of the filtered predictors along the
block's coordinates, columns of the forward Jacobian of the recursion.
The level of the subformula is its unpenalized sub-terms, read from
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html),
and the constant direction is the combination \\\gamma\\ of the level's
coordinates with \\Z\_{\rm level}\gamma = 1\\, the derivative of the
predictors along a constant shift of the parameter. The span is tested
on the subformula's own design, by
[`jfp_constant()`](https://statmodels7.github.io/statmodels7/reference/jfp_constant.md).
