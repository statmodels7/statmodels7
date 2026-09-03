# One Evaluation Point, Shared

The quantities the marginal criterion, its gradient and its Hessian all
read at one \\(\beta, \theta)\\: the linear predictors, the information,
the penalty's Hessian, their sum and its factorization.

## Usage

``` r
outer_context(spec, design, coef, hyper, approx = "opg")
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- approx:

  The approximation for the expected information.

## Value

An environment carrying the point and, as they are asked for, the
quantities derived from it.

## Details

The three consumers each used to assemble these for themselves, so at
one point the information was built three times and the penalized matrix
factorized twice, and
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
additionally calls the gradient, which repeated the whole of it a fourth
time. Measured by `Rprof`'s `by.total` on a random intercept over 500
levels at 20000 observations, the gradient and the Hessian together
accounted for 128 per cent of the fit, the overlap being exactly that
repetition.

The context is an environment, so an accessor fills it in place and
later readers find the quantity already there. Each is computed on first
demand and never speculatively: the criterion alone does not need an
inverse, and a search running without an exact gradient must not pay for
one.

Passing `NULL` wherever a context is accepted restores the earlier
behavior exactly, so every existing caller goes on working unchanged, a
caller's own code included.

## See also

[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`statmod_marginal_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_hess.md)
