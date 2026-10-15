# Which Optimizer the Outer Search Uses When the Caller Names None

The choice is made from what the criterion can supply: its exact
Hessian, its exact gradient, or neither.

A covariance class split between the coefficients and a filter's own
parameters takes
[`optimizers7::lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.html)
even where the exact gradient is unavailable, because such a class
carries at least three hyperparameters (two scales and a correlation)
and a simplex search is unreliable at that dimension.

[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
is given `max_length = 5`, which bounds every component of a step on the
free scale, as `maxNstep` does in mgcv. Where the criterion flattens as
a smoothing parameter grows, the curvature at the start can be small
enough for the Newton step to jump past the maximum onto that plateau,
where the gradient is close to zero and the search stops.

A search over one coordinate with no exact gradient takes
[`optimizers7::brent()`](https://statmodels7.github.io/optimizers7/reference/brent.html),
which brackets the maximum and reduces the bracket by parabolic steps,
stopping when its half-width on the free scale is below \\10^{-3}\\, a
change of one part in a thousand in the hyperparameter. On a single
coordinate it needs fewer evaluations than
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html).

## Usage

``` r
outer_default_optimizer(exact, use_hess, mixed = FALSE, n = NA_integer_)
```

## Arguments

- exact:

  Whether the criterion has an exact gradient.

- use_hess:

  Whether the search should steer by the exact Hessian, which is not the
  same question as whether one exists: see
  [`outer_newton_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_newton_ok.md).

- mixed:

  Whether the model carries a covariance class spanning the coefficients
  and a structural term's own parameters.

- n:

  The number of coordinates the search moves.

## Value

An optimizers7 optimizer.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
