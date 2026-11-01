# The Movement of the Third and Fourth Derivatives of a Jeffreys Block

`jfp_mv_ctx()` sets up, for one block and up to two directions, the
context from which `jfp_mv_white()`, `jfp_mv_B()` and `jfp_mv_pair()`
compute the three objects of the third and fourth derivatives of
\\\log\|\tilde I\|\\: the whitened derivative of the block's matrix
along the directions, the rows \\B_R\\ of the low-rank terms, and the
trace terms \\T_R(Y)\_{kl} = \partial_R\\\mathrm{tr}(\tilde
Y\\\partial\_{kl} \tilde I)\\ with \\\tilde Y = U^{-1}YU^{-\top}\\.

## Usage

``` r
jfp_mv_ctx(st, bs, ds)

jfp_mv_white(cx, R)

jfp_mv_B(cx, R)

jfp_mv_pair(cx, Y, R)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- bs:

  A block's entry of the state.

- ds:

  A named list of directions (`v`, and `w` for the fourth derivative),
  each from
  [`jfp_dir()`](https://statmodels7.github.io/statmodels7/reference/jfp_dir.md).

- cx:

  A context, from `jfp_mv_ctx()`.

- R:

  A character vector naming the directions taken, a subset of
  `names(ds)`.

- Y:

  A symmetric \\\tilde p\times\tilde p\\ matrix.

## Value

`jfp_mv_ctx()` an environment; `jfp_mv_white()` a \\\tilde p\times\tilde
p\\ matrix; `jfp_mv_B()` a matrix over the joint vector and the \\\tilde
p^2\\ entries; `jfp_mv_pair()` a list with `sparse` (per-observation
weights keyed as
[`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
reads its `H`) and `dense` (a list of entries `idx`, `X`, each adding
`X` in the rows `idx` and its transpose in the columns `idx`).

## Details

With the block's rows in one equation, \\\tilde I = \sum_t w_tJ_t \tilde
V_t^\top\tilde V_t\\ with \\J_t\\ the family's expected information of
that equation. Each derivative distributes the directions over \\J\\ and
the two factors \\\tilde V\\ by the Leibniz rule
([`jfp_leibniz()`](https://statmodels7.github.io/statmodels7/reference/jfp_leibniz.md)),
and a derivative of \\J\\ is the Faa di Bruno sum in the predictors,
with the second derivative of a predictor where a term in its equation
recomputes its block at the coefficients
([`jfp_mv_J()`](https://statmodels7.github.io/statmodels7/reference/jfp_dir.md)).
A free index enters the information through the Jacobian column and
through the derivative of the moving block along the directions, and it
enters the rows only in the term that holds the block. For a model whose
predictors are linear and a block whose rows do not move, every share
beyond the design's cross products is zero and the objects are those of
the static derivatives.

In a model with a filter the Jacobian columns are the joint rows of
every equation
([`jfp_jac()`](https://statmodels7.github.io/statmodels7/reference/jfp_dir.md)),
and each equation that the filter drives is a moving unit whose
derivatives along the directions come from the filter's kernels
([`jfp_mv_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_dir.md)).
A block in the subformula of a filter's parameter has rows in every
equation that the filter drives, so \\\tilde I = \sum_t
w_t\sum\_{q,r}I\_{qr,t}\tilde V\_{q,t}^\top \tilde V\_{r,t}\\ and the
three objects are summed over the ordered pairs \\(q, r)\\
([`jfp_mv_setpair()`](https://statmodels7.github.io/statmodels7/reference/jfp_dir.md));
the rows' derivatives are read along the block's whitened directions
\\d_i\\, one kernel call per direction.
