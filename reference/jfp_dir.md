# The Building Blocks of the Moving Shares

`jfp_dir()` carries a direction over the joint vector with the movement
of every predictor along it. `jfp_mv_J()` returns the derivative of the
expected information of the block's equation, with the predictor indices
`fixed`, along the directions `R`. `jfp_mv_e2()` returns the second
derivative of every predictor along the two directions. `jfp_mv_W()` and
`jfp_mv_Wk()` return the whitened derivatives of the block's rows along
`R`, the second one also in each coefficient of the term that holds the
block. `jfp_mv_etak()` returns the derivative of a moving term's block
along `R`; `jfp_mv_etakk()` its derivative in two of the term's
coefficients along `R`, contracted against per-observation weights;
`jfp_mv_cv()` the derivative of the block's rows in two coefficients of
its term along `R`, contracted against weights on the block's columns.
`jfp_mv_pfun()` wraps a product of two factors with a memory, and
`jfp_nonempty()` lists the non-empty subsets of at most two directions.
In a model with a filter, `jfp_mv_setpair()` selects the ordered pair of
equations whose information and rows the other functions read,
`jfp_mv_full()` places directions over the filter's full width,
`jfp_mv_rows()` returns the predictors' derivative along the directions
(and along the block's whitened direction `i`, when given), one row per
observation and per equation, and `jfp_jac()` returns an equation's
Jacobian columns with their positions in the joint vector; the weights
`Q` of `jfp_mv_cv()` are then on the whitened directions.

## Usage

``` r
jfp_dir(st, v)

jfp_mv_J(cx, fixed, R)

jfp_mv_e2(cx)

jfp_mv_W(cx, S, eq = cx$a)

jfp_mv_Wk(cx, S, eq = cx$a)

jfp_mv_etak(cx, i, S)

jfp_mv_etakk(cx, i, S, g)

jfp_mv_cv(cx, S, Q, eq = cx$a)

jfp_mv_pfun(pr, A, B)

jfp_nonempty(R)

jfp_mv_setpair(cx, k)

jfp_mv_full(cx, S)

jfp_mv_rows(cx, i, S)

jfp_jac(st, cc)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- v:

  A direction over the joint vector.

- cx:

  A context, from
  [`jfp_mv_ctx()`](https://statmodels7.github.io/statmodels7/reference/jfp_mv_ctx.md).

- fixed:

  Equation indices of the derivative in the predictors.

- R, S:

  A character vector naming directions.

- eq:

  The equation of a factor of the block's rows.

- i:

  The position of a moving term among the context's units.

- g:

  Per-observation weights.

- Q:

  Weights over the observations and the block's columns.

- pr:

  A product of two matrices.

- A, B:

  Functions of a set of directions returning a factor.

- k:

  The position of an ordered pair of equations in the context.

- cc:

  An equation.

## Value

As described; `NULL` where a factor is zero.
