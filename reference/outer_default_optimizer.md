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

## Usage

``` r
outer_default_optimizer(exact, use_hess, mixed = FALSE)
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

## Value

An optimizers7 optimizer.

## See also

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
[`outer_gradient_ok()`](https://statmodels7.github.io/statmodels7/reference/outer_gradient_ok.md)
