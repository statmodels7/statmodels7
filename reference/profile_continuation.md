# Profile by Continuation

Selects a method of
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
that follows the refitted coefficients from one value of the grid to the
next. At each value a predictor extrapolates the coefficients that are
not held from the previous point, and a Newton corrector brings them to
the optimum with the held coefficients fixed.

## Usage

``` r
ProfileContinuation(
  name = character(0),
  order = integer(0),
  maxit = integer(0)
)

profile_continuation(order = 1, maxit = Inf)
```

## Arguments

- name:

  The short name of the method, a single string. The class constructor
  `ProfileContinuation()` takes it, and `profile_continuation()` sets it
  to `"continuation"`.

- order:

  The order of the predictor: 1 for the tangent, 2 for the tangent and
  the curvature.

- maxit:

  The number of passes the corrector may take at each point, or `Inf`
  for convergence.

## Value

A `ProfileContinuation` object, inheriting from
[profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md),
with the properties `name`, `order` and `maxit`.

## Details

The predictor of order 1 is the tangent of the trace of the refitted
coefficients, \\d\hat\eta/d\theta = -H\_{\eta\eta}^{-1}
H\_{\eta\theta}\\, where \\\theta\\ are the held coefficients, \\\eta\\
the others and \\H\\ the Hessian of the objective of the fit (the
negative penalized log-likelihood) at the previous point. Order 2 adds
the curvature of the trace, obtained from the second derivative along
the step of the gradient of the objective with the trace developed to
order 1, by one central stencil. The corrector is the refit of
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md)
started from the prediction. With `maxit = Inf` it runs to convergence,
and the rows equal those of
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md)
within the tolerance of the fit. A finite `maxit` limits the passes of
the corrector, and the rows are then approximations whose `converged` is
`FALSE` where the limit stopped the corrector.

If the corrector signals an error, returns a non-finite objective or,
with `maxit = Inf`, does not converge, the step is halved: a point
halfway between the previous point and the requested one is computed
first, and the requested point is attempted again from it. The
intermediate points are not part of the output.

Coordinates that the outer criterion estimates (the unpenalized
coefficients of the parameters other than the first under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md))
are re-estimated by the criterion at every point of the exact method.
They are followed through the nested problem: the held coefficients and
those coordinates are criterion coordinates together, the exact Hessian
of the criterion over them gives the derivative of the estimated
coordinates with respect to the held ones, and the other coefficients
are the penalized mode at the predicted coordinates. The first step
starts from the refit with the held coefficients at the estimate, which
differs from the fit where the criterion does not integrate the held
coefficients.

The method requires `hyper = "held"` and a model whose coefficients are
smooth functions of the held ones: the terms
[`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html),
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html),
[`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html),
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html),
[`modelterms7::ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
and the smoothed forms of
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html).
A model with a break-point term that is not smoothed, a kinked penalty
or a structural term is rejected with a pointer to
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md).

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md),
[`profile_taylor()`](https://statmodels7.github.io/statmodels7/reference/profile_taylor.md).

## Examples

``` r
profile_continuation()
#> <statmodels7::ProfileContinuation>
#>  @ name : chr "continuation"
#>  @ order: int 1
#>  @ maxit: num Inf
profile_continuation(order = 2, maxit = 1)
#> <statmodels7::ProfileContinuation>
#>  @ name : chr "continuation"
#>  @ order: int 2
#>  @ maxit: num 1
```
