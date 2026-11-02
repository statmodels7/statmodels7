# A Recomputed Block Beside a Structural Term, in the Outer Hessian

Returns what the exact outer Hessian of the structural branch
([`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md))
needs, beyond
[`structural_refresh_move()`](https://statmodels7.github.io/statmodels7/reference/structural_refresh_move.md),
where a term recomputes its block at the coefficients: the movement of
\\K\\ with the block along each hyperparameter's movement of the mode,
the second movement of \\K\\ along two of them, and the movement of the
block's term of the likelihood's curvature.

## Usage

``` r
structural_refresh_second(spec, design, coef, jd, M, rm, bhat)
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

- M:

  The matrix against which the traces are taken.

- rm:

  The result of
  [`structural_refresh_move()`](https://statmodels7.github.io/statmodels7/reference/structural_refresh_move.md).

- bhat:

  The mode's movement along each hyperparameter, over the estimated
  coordinates.

## Value

`NULL` where the wider information cannot be formed; otherwise a list
with `R`, the matrices \\\partial_R K\[D_m\]\\; `mode(m, l)`, the vector
\\\partial D^{\rm m}\[b_l\]\\b_m\\; and `trace(m, l)`, the trace against
\\M\\ of the second movement of \\K\\ that the rows carry.

## Details

Write \\K = K(u, R)\\ for the information at the joint vector \\u\\ with
static rows \\R\\, which the block makes a function of \\u\\: \\\partial
R/\partial u\\\[b\] = D_b\\
([`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)),
and the static predictor's second derivative is \\D_b c\\. At fixed
predictors and filter parameters \\K\\ is a quadratic form in the rows,
\\K\_{ab} = B(R_a, R_b)\\, where a column of the filter's own parameters
enters linearly. Adding the columns of every \\D_m\\ to the rows as
further unknowns, placed before the filter's parameters, gives an
information \\K^+\\ whose blocks are the values of \\B\\ on those
columns, so that \$\$\partial_R K\[D\] = \sum_j (e_j k_j^\top + k_j
e_j^\top),\qquad \partial^2_R K\[D_m, D_l\]\_{ab} = B(D\_{m,a},
D\_{l,b}) + B(D\_{l,a}, D\_{m,b}),\$\$ with \\k_j\\ the column of
\\K^+\\ for the added unknown that carries column \\j\\ of \\D\\. The
derivative of \\K^+\\ along \\b\\ at fixed rows
([`structural_dk_matrix()`](https://statmodels7.github.io/statmodels7/reference/structural_dk_matrix.md)
on the wider rows) gives the mixed derivative \\\partial_u\partial_R
K\[b_m, D_l\]\\, and its derivative along an added unknown gives the
movement of \\K\\ with the static predictor's second derivative \\D_m
b_l\\. The term \\\partial_R K\[\partial^2 R\[b_m, b_l\]\]\\, with the
rows' second derivative from
[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html),
is the exact central difference of the quadratic.

The block's term of the curvature, \\(D^{\rm m} b)\_c = -\sum_t g_t\\
(D_b)\_{t,c}\\ with \\g\\ the total derivative in the static predictor
([`statmod_score_obs()`](https://statmodels7.github.io/statmodels7/reference/statmod_score_obs.md)),
moves along \\b_l\\ through \\g\\, which is the row of \\K^+\\ for the
added unknown, and through the rows' second derivative.

## See also

[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md),
[`structural_refresh_move()`](https://statmodels7.github.io/statmodels7/reference/structural_refresh_move.md)
