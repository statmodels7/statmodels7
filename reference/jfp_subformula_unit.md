# A Jeffreys Block in the Subformula of a Refreshed Term

Resolves a Jeffreys unit that a term such as
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
or a smoothed
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
carries up from the subformula of one of its parameters: the term's
columns from which the block's rows are combined, and whether the
constant direction of the subformula enters them.

## Usage

``` r
jfp_subformula_unit(u, spec, design)
```

## Arguments

- u:

  The unit.

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

## Value

A list with `kind` (`"refresh"`), `const`, `basis` and `C`, as
[`jfp_block()`](https://statmodels7.github.io/statmodels7/reference/jfp_block.md)
takes them.

## Details

The level of the subformula is its unpenalized sub-terms, read from
[`modelterms7::term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.html).
The constant direction is the derivative of the predictor along a
constant shift of the parameter, which is the combination \\X\_{\rm
level}\gamma\\ of the level's refreshed columns with \\Z\_{\rm
level}\gamma = 1\\; the span is tested on the subformula's own design,
by
[`jfp_constant()`](https://statmodels7.github.io/statmodels7/reference/jfp_constant.md).

A term whose block is a working linearization (a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html))
is rejected, and so is the break-point subformula of a sharp
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
whose rows jump when the break-point crosses an observation. A parameter
whose own development moves is rejected as well.
