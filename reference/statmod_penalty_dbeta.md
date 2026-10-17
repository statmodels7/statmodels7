# How a Penalty's Hessian Moves With the Mode

Computes \\\sum_u \partial S_u/\partial\beta\\\[v\]\\ over every penalty
on the stacked coefficients whose Hessian depends on them, placed in the
stacked coefficient space. `penalty_dbeta_trace()` computes its trace
against a matrix, and `penalty_dbeta_blocks()` returns the per-penalty
pieces from which both are assembled.

## Usage

``` r
statmod_penalty_dbeta(spec, design, coef, hyper, v, total)

penalty_dbeta_trace(spec, design, coef, hyper, v, M)

penalty_dbeta_blocks(spec, design, coef, hyper, v)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- v:

  The direction, over the stacked coefficients.

- total:

  The stacked width.

- M:

  The matrix against which the trace is taken, `total` square.

## Value

`statmod_penalty_dbeta()` a `total` square matrix, or `NULL` if no
penalty's Hessian moves; `penalty_dbeta_trace()` a single number;
`penalty_dbeta_blocks()` a list of `pos` and `T` pairs.

## Details

A marginal criterion's determinant is of \\H + S\\, so the mode's
movement \\v\\ reaches it through \\S\\ wherever \\S\\ depends on the
coefficients, as it does for a heavy-tailed prior on a random effect,
and the piece is \\\mathrm{tr}(M\\\partial S/\partial\beta\[v\])\\. A
prediction-error criterion reads the matrix itself inside its trace.
Every such penalty contributes, not only the one owning the
hyperparameter being differentiated, because the mode moves in all of
the coefficients at once.

A penalty for which
[`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html)
returns `TRUE` is skipped (this is every ridge, smooth and Gaussian
random effect), so a model carrying only those assembles nothing: the
matrix is `NULL` and the trace exactly `0`, which leaves the gradient
unchanged. A kinked penalty and a mixed covariance class are skipped as
well, and a penalty over a structural term's own parameters is read on
the joint vector by
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md).

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_pe_derivs()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe_derivs.md),
[`penalties7::penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.html)
