# The Joint Second Derivative of the Penalized Information, Traced

\\\mathrm{tr}(M\\\partial^2 K/\partial u^2\[v, w\])\\: the quantity
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
gives at the first order, one order up and contracted against a second
direction.

## Usage

``` r
structural_chain_extra2(spec, design, jd, M, st, blk4, v, w)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows.

- M:

  The matrix the trace is taken against.

- st:

  The shared quantities, from
  [`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md).

- blk4:

  The `blocks` factory carrying the family's fifth derivative, as
  [`.structural_blocks()`](https://statmodels7.github.io/statmodels7/reference/dot-structural_blocks.md)
  builds it.

- v, w:

  The two directions, over the estimated coordinates.

## Value

A single number.

## Details

Nine terms, and none of them assembles a matrix over the coefficients: a
term in \\V_a^\top V_b\\ traces as a weighted sum of the per-observation
diagonal \\G\\, a term carrying \\\mathrm{d}\varphi\\ traces against the
rows \\M\\ has already been applied to, and the three terms carrying the
recursion's own derivatives trace against what
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html),
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
and
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html)
return at the right weights.

The last of those is the fourth derivative of the predictor through the
recursion, which is the object this whole route exists for, and the only
place the family's FIFTH derivative enters is `P` inside its `blocks`
callback.

## See also

[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md),
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
