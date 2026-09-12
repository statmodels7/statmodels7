# The Outer Criterion's Curvature at a Reported Point

\\A = -(H + H^\top)/2\\, the matrix a Newton decrement is read against,
from the analytic Hessian where the form has one and from one central
difference of the **exact** gradient where it does not.

## Usage

``` r
outer_curvature(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  inner = NULL
)
```

## Arguments

- spec, design, coef, hyper, method, idx, basis:

  As
  [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
  takes them.

- inner:

  The inner optimizer, which the differenced route refits with.

## Value

A list with `A`, the symmetric negated Hessian or `NULL`; `source`,
`"analytic"` or `"differenced"`; and `why`, the reason where there is no
`A`.

## Details

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads the rise the criterion would still buy, which needs a curvature
and not only a gradient. The analytic route is
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
gated by
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
at order 2 – the gate is load-bearing rather than defensive, since the
assembly returns a number wherever it is called and the wrong one where
the order is not covered.

**Where the order is not covered the gradient is differenced instead**,
which is
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
and the same route
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
already takes for a model carrying a filter. The two forms that reach it
are a criterion asked for on the **expected** information, where order 2
is refused because the criterion's own second derivative would want the
next order of \\\partial\mathbb{E}\[\ell''\]/ \partial\eta\\, and a
**separable** penalty whose curvature moves with the coefficient, a
heavy-tailed prior on a random effect among them. ⚠️ The second is not a
missing derivative: measured, a t prior answers
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html),
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html)
and
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html).
What refuses it is
[`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html),
TRUE for a ridge and a gaussian prior and FALSE here, the order-2
assembly being written for a penalty whose Hessian in the coefficients
does not move with them. Both forms carry an exact gradient, which is
what the difference is taken of.

The two routes agree where both exist: on a gaussian smooth the
decrement reads `3.01e-10` by either. What the difference costs is
measured and is not nothing – it is \\4n_h\\ refits, 0.05 to 0.13
seconds on the three shapes tried, 5 to 37 per cent of the fit itself,
where the analytic route is below the clock's own resolution. So a fit
whose form carries the analytic Hessian refits nothing, and one that
does not pays for its verdict once, at the summary.

⚠️ The difference is a stopgap and not the destination. One gap is a
quantity nobody has written – the second derivative of an expected
information – and the other is an assembly written under an assumption
one penalty does not satisfy; until both are closed, a certificate on
those forms rests on a stencil where every other reading in this package
rests on algebra.

⚠️ `source` is informative and is not load-bearing, and one case can
mislabel: for a model carrying a filter
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
falls back to
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
itself where the analytic assembly fails, and that is reported here as
`"analytic"`, this function having asked only whether the order was
covered. Nothing reads `source` to decide anything.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
