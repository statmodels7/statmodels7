# The Outer Criterion's Curvature at a Reported Point

Computes \\A = -(H + H^\top)/2\\, the matrix against which a Newton
decrement is read, from the analytic Hessian where the form has one and
from one central difference of the exact gradient where it does not.

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
  inner = NULL,
  gam = NULL
)
```

## Arguments

- spec, design, coef, hyper, method, idx, basis:

  As
  [`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
  takes them.

- inner:

  The inner optimizer with which the differenced route refits.

- gam:

  The coefficients that the criterion estimates, from
  [`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
  or `NULL`.

## Value

A list with `A`, the symmetric negated Hessian or `NULL`; `source`,
`"analytic"` or `"differenced"`, or `NA` where there is no `A`; and
`why`, the reason where there is no `A`.

## Details

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
reads the rise that the criterion would still buy, which needs a
curvature and not only a gradient. The analytic route is
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
gated by
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
at order 2. The gate is required for correctness, since the assembly
returns a number wherever it is called and a wrong one where the order
is not covered.

Where the order is not covered, the exact gradient is differenced
instead by
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md),
which is also the route that
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
takes for a model carrying a filter. The forms that reach it are a
criterion on the expected information for a family that does not supply
the second derivative of its expected information (or for a model with a
refreshable block or a structural term), and a penalty whose Hessian
moves with the coefficients, such as a heavy-tailed prior on a random
effect, placed beside a structural term, where
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
reads no movement of that Hessian. The same prior without a structural
term is analytic, through
[`statmod_penalty_second()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_second.md).
Both kinds of form carry an exact gradient, and the difference is taken
of that gradient.

The two routes agree where both exist. The difference costs \\4(n_h +
n\_\gamma)\\ refits, where \\n\_\gamma\\ is the number of coefficients
that the criterion estimates; the analytic route refits nothing.

`source` is informative only. For a model carrying a filter
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
falls back to
[`statmod_hess_stencil()`](https://statmodels7.github.io/statmodels7/reference/statmod_hess_stencil.md)
where the analytic assembly fails, and `source` then still reads
`"analytic"`, because this function has checked only whether the order
was covered.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`joint_decrement()`](https://statmodels7.github.io/statmodels7/reference/joint_decrement.md)
