# The Third and Fourth Derivatives of the Jeffreys Penalty

`jfp_deriv3_dir()` returns the derivative of the Hessian of \\\sum_b
a_b\log\|\tilde I_b\|\\ along a direction \\v\\ of the joint vector, and
`jfp_deriv4_dir()` its second derivative along two directions \\v\\ and
\\w\\, as matrices over the joint vector. `jfp_deriv3()` and
`jfp_deriv4()` return their traces against a symmetric matrix \\M\\,
without assembling them.

## Usage

``` r
jfp_deriv3_dir(st, v)

jfp_deriv3(st, M, v)

jfp_deriv4_dir(st, v, w)

jfp_deriv4(st, M, v, w)
```

## Arguments

- st:

  A state, from
  [`jfp_state()`](https://statmodels7.github.io/statmodels7/reference/jfp_state.md).

- v, w:

  Directions over the joint vector.

- M:

  A symmetric matrix over the joint vector.

## Value

`jfp_deriv3_dir()` and `jfp_deriv4_dir()` a base matrix over the joint
vector; `jfp_deriv3()` and `jfp_deriv4()` a single number.

## Details

The derivatives are those of the determinant with the derivatives of the
family's expected information to the fourth order, from
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and its siblings of orders two to four. A term with both indices inside
one derivative of the block's matrix is a weighted cross product of the
designs and has the sparsity of the information; the others are products
of low rank through matrices of \\\tilde p^2\\ columns. Where a term
recomputes its block at the coefficients
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
a smoothed
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)),
the derivatives of the predictors to the fourth order and, for a block
in the subformula of such a term, the movement of the block's rows enter
as well
([`jfp_mv_ctx()`](https://statmodels7.github.io/statmodels7/reference/jfp_mv_ctx.md));
the derivatives of the term's block come from
[`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)
and its siblings up to
[`modelterms7::term_block_deriv4()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv4.html).
In a model with a filter the predictors are those of the joint vector
([`jfp_joint_rows()`](https://statmodels7.github.io/statmodels7/reference/jfp_joint_rows.md)):
their derivatives come from the filter's kernels along the directions,
and for a block in the subformula of a filter's parameter the rows'
derivatives are read along the whitened directions of the block, one
kernel call per direction
([`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html),
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html),
[`modelterms7::term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.html)).

The outer criterion reads them: the derivative of its penalized matrix
along the movement of the mode carries \\-a\\D_3(v)\\, and its second
derivative \\-a\\D_4(v, w)\\
([`statmod_penalty_dbeta()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_dbeta.md),
[`statmod_penalty_second()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_second.md)).

## See also

[`jfp_hessian()`](https://statmodels7.github.io/statmodels7/reference/jfp_hessian.md)
