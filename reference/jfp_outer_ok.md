# Whether the Outer Derivatives Can Read a Jeffreys Penalty

Returns `TRUE` where the exact outer gradient (`order = 1`) or Hessian
(`order = 2`) can carry the movement of a Jeffreys penalty.

## Usage

``` r
jfp_filter_outer_ok(spec, design, order = 1L)

jfp_outer_ok(spec, design, order = 1L)

jfp_refresh_ok(spec, design, order = 1L)

jfp_outer_probe(d, order)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- order:

  `1` or `2`.

- d:

  A distribution.

## Value

A single logical.

## Details

The gradient reads the third derivative of the penalty, so the third
derivative of the family's expected information, and the Hessian the
fourth; each is probed once at a random point with the method that
[`jfp_emethod()`](https://statmodels7.github.io/statmodels7/reference/jfp_emethod.md)
selects.

Beside a filter
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
a single or a joint one) the derivatives are written on the joint
vector, which is the vector of the criterion only where a penalty covers
the filter's own parameters
([`structural_penalized()`](https://statmodels7.github.io/statmodels7/reference/structural_penalized.md));
elsewhere the function returns `FALSE`. The gradient then reads the
filter's third derivative and the family's fourth
([`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)),
the Hessian one order more of each; a block in the subformula of a
filter's parameter reads one order more again, the Hessian through
[`modelterms7::term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.html),
with the fifth derivative of the expected information where the score is
scaled
([`distributions7::distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.html)).
A filter together with a term whose block is recomputed at the
coefficients is not covered.

Where a term recomputes its block at the coefficients,
`jfp_refresh_ok()` asks that the block be a Jacobian
([`modelterms7::term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.html))
and that the term's class have its own method for the derivative of its
block that the order reads:
[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html)
for the gradient and
[`modelterms7::term_block_deriv3()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv3.html)
for the Hessian, one order more for the term that holds a Jeffreys block
in its subformula, whose rows move. The method is called once at the
term's coefficients, so that a term that signals an error there (an
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
parameter whose development moves, a step's fourth derivative) is read
as `FALSE`.
