# A Recomputed Block Beside a Structural Term, in the Outer Gradient

Returns what the exact outer gradient of the structural branch
([`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md))
needs where a term recomputes its block at the coefficients
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
a smoothed
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)):
the block's second-derivative term of the likelihood's curvature, by
which the mode moves, and the movement of \\K\\ with the block.

## Usage

``` r
structural_refresh_move(spec, design, coef, jd)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design, refreshed at `coef`.

- coef:

  The coefficients at the mode.

- jd:

  The joint rows, from
  [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md).

## Value

`NULL` where no term recomputes its block; otherwise a list with `Dm`, a
matrix over the coefficients, and `trace`, a function of the direction
over the joint vector and of \\M\\.

## Details

The criterion's \\K\\ is assembled on the block at the mode and, as on
the route without a structural term, leaves out the term \\\sum_t\sum_q
g\_{q,t}\\\partial^2\eta^{\rm s}\_{q,t}/\partial u^2\\ that the second
derivative of the static predictors \\\eta^{\rm s}\\ contributes to the
likelihood's curvature. The mode moves by the full curvature, so `Dm` is
that term, from
[`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md)
with \\g\_{q,t}\\ the total derivative of the log-likelihood in the
static predictor (the score plus the filter's reverse pass,
[`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md));
the filter carries the block's second derivative into later levels
through that reverse pass.

Along a direction \\v\\ the static rows move by \\D = \partial
X/\partial u\\\[v\]\\
([`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)),
and \\K\\ moves through the filter's Jacobian and curvature, which the
recursion computes from those rows. At fixed predictors and filter
parameters \\K\\ is a polynomial of degree two in the rows: the Jacobian
is affine in them, and the curvature and \\\sum H V^\top V\\ are
quadratic. Its derivative along \\D\\ is therefore exactly \\(K(V +
sD) - K(V - sD))/(2s)\\ for any \\s\\, which
[`trace()`](https://rdrr.io/r/base/trace.html) takes against \\M\\ with
\\s\\ scaling \\D\\ to the size of \\V\\. The movement of \\K\\ through
the predictors and the filter's parameters is read by
[`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md)
and
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md).

## See also

[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
the same movement on the route without a structural term.
