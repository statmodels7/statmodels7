# The Outer Hessian by One Difference of the Exact Gradient

The criterion's second derivative in the hyperparameters where no
analytic route exists: one central difference of
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
with the coefficients refitted at every probe.

## Usage

``` r
statmod_hess_stencil(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  inner = NULL,
  h = hess_stencil_step()
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients the probes start from.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The hyperparameter index, from
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- basis:

  The integrated basis, or `NULL`.

- inner:

  The inner optimizer the probes refit with;
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  where none is given.

- h:

  The step, on the free scale the search runs on.

## Value

The Hessian on the free scale, or `NULL` where a probe could not be
evaluated or the two steps disagree.

## Details

This is what a model carrying a structural term gets, and the reason is
that the analytic assembly of
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
spans the stacked coefficients while a filter's own parameters are
estimated beside them and move with the hyperparameter as well.
Extending it is not a matter of bookkeeping: each order of
differentiation through the recursion pulls in one more order of the
response's family, so the first derivative reads the family's fourth
through
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
and the second would read a fifth, which does not exist.

Differencing an ANALYTIC quantity once is the licence this toolkit
already grants itself for the Student t's degrees of freedom and for the
marginal break-point's prior rows. What it forbids is a difference of a
difference, and there is one layer here.

## The two probes start from the same place

Both refit from the coefficients given and from the structural state as
it stands, restored before each probe, so the two differ in the
hyperparameter alone. That is what makes the result stable: the mode's
own location error is nearly the same at \\+h\\ and \\-h\\ and cancels
in the difference rather than being amplified by \\1/h\\. Measured on a
penalized filter against a second difference of the criterion, the
result is FLAT over four decades of the step – -1.971435 at every \\h\\
from 1e-2 down to 3e-5, against a criterion second difference of
-1.971456 – where the assembled Hessian reads -1.1e-08.

## Where it refuses

A stencil inherits the reproducibility of the quantity it differences,
divided by the step, so where the mode is poorly located the answer is
noise rather than a curvature. It is therefore computed TWICE, at `h`
and at `3 * h`, and refused where the two disagree by more than
[`hess_stencil_tol()`](https://statmodels7.github.io/statmodels7/reference/hess_stencil_step.md).
The two regimes are six orders apart and nothing sits between them:
measured, 5.5e-08 on a penalized filter, 9.3e-07 on an unpenalized one
beside a smooth and 1.1e-06 on an ordinary smooth, against 6.0e-01 and
5.7e-01 on two mixed covariance classes whose correlation the search
left at \\\|z\| = 8.53\\, where the chart's conditioning is 1e10 and the
exact gradient itself reads 1e-3. A refusal there is the answer, and the
consumers already handle `NULL`.

## See also

[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`hess_stencil_step()`](https://statmodels7.github.io/statmodels7/reference/hess_stencil_step.md)
