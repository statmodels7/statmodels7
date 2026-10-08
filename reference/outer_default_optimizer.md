# Which Optimizer the Outer Search Uses When the Caller Names None

The choice is made from what the criterion can supply: its exact
Hessian, its exact gradient, or neither.

A covariance class split between the coefficients and a filter's own
parameters takes `lbfgs()` even where the exact gradient is unavailable
to it, which since the gradient was written is a shape it reaches only
for some further reason of its own. Such a class carries at least three
hyperparameters – two scales and a correlation – and three is where this
file already records the simplex stalling. Measured on a panel of twelve
groups while the gradient was still refused, both searches reporting
convergence: `nelder_mead()` reached a criterion of -684.789 in 123
seconds, at a loading scale collapsed to 2.1e-06 and a correlation of
-0.9999, while `lbfgs()` differencing the criterion reached -683.305 in
27 seconds at scales of 0.647 and 0.389 and a correlation of 0.733 –
which on data simulated with the two effects INDEPENDENT is the answer
the two separate fits give, to 0.2 of log-likelihood.

`newton()` is given `max_length = 5`, which bounds every component of a
step on the free scale, as `maxNstep` does in mgcv. Where the criterion
flattens as a smoothing parameter grows, the curvature at the start can
be small enough for the Newton step to jump past the maximum onto that
plateau, where the gradient is close to zero and the search stops.
Measured on a ridge over the fifteen predictors of
[`MASS::UScrime`](https://rdrr.io/pkg/MASS/man/UScrime.html),
standardized by the term, the first step went from \\\lambda = 1\\ to
\\1.6\times10^{10}\\, a length of 23.5, and the fit reported the empty
model at a criterion of -65.40, where the maximum is -53.30 near
\\\lambda = 91.5\\.

A search over one coordinate with no exact gradient takes
[`optimizers7::brent()`](https://statmodels7.github.io/optimizers7/reference/brent.html),
which brackets the maximum and reduces the bracket by parabolic steps,
stopping when its half-width on the free scale is below \\10^{-3}\\, a
change of one part in a thousand in the hyperparameter. Measured on
`seg(t, psi ~ random(~1 | id), marginal = TRUE)` over a 20 x 15 panel,
where the only hyperparameter is the scale of a random intercept, the
search took 11 evaluations against the 23 of
[`optimizers7::nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.html),
at a criterion 2e-05 lower.

## Usage

``` r
outer_default_optimizer(exact, use_hess, mixed = FALSE, n = NA_integer_)
```

## Arguments

- exact:

  Whether the criterion has an exact gradient.

- use_hess:

  Whether the search should STEER by the exact Hessian, which is not the
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
