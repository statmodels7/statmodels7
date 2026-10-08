# The Quantities of a Structural Term That Read a Coordinate at an Edge

Marks each quantity
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)
reports for a structural term according to whether it depends on a free
coordinate that has run past the edge of a chart mapping onto a bounded
set.

## Usage

``` r
structural_edge_rows(tm, zeta, nm, edge = 8)
```

## Arguments

- tm:

  The structural term.

- zeta:

  Its parameters on the unconstrained scale.

- nm:

  Its parameter names, as
  [`modelterms7::term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.html)
  gives them.

- edge:

  The free value past which a coordinate is at the edge.

## Value

A logical vector with one element per reported quantity.

## Details

The coordinates checked are those
[`modelterms7::term_charted()`](https://statmodels7.github.io/modelterms7/reference/term_charted.html)
names, and a coordinate is at the edge when its free value exceeds
`edge` in absolute value, the rule
[`edge_violations()`](https://statmodels7.github.io/statmodels7/reference/edge_violations.md)
and
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
use, on a side where its chart maps onto a finite bound: a log link has
an edge towards zero and none towards infinity. The dependence is read
off the Jacobian evaluated with those coordinates moved back to
\\\pm\\`edge`: at the fitted point the derivative of a saturated chart
can underflow to zero, which would hide the dependence it is meant to
reveal.

## See also

[`statmod_structural_table()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_table.md),
[`readable_joint()`](https://statmodels7.github.io/statmodels7/reference/readable_joint.md)
