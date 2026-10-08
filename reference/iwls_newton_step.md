# One Newton Step on the Full Hessian and Its Line Search

Solves for the Newton increment on a given Hessian over the free
coordinates and halves the step until Armijo's condition holds or the
budget of halvings is spent.

## Usage

``` r
iwls_newton_step(obj, beta, value, g, H, frozen, halvings)
```

## Arguments

- obj:

  The objective.

- beta:

  The current coefficients.

- value:

  The objective there.

- g:

  The gradient there.

- H:

  The full Hessian there.

- frozen:

  The held positions.

- halvings:

  The budget of step halvings.

## Value

`NULL` where the Hessian is not usable, and otherwise a list of `ok`,
`cand`, `vnew` and `step_used`.

## Details

The Hessian of an objective whose design moves with its coefficients
need not be positive definite away from the mode, so its eigenvalues are
replaced by their absolute values and floored at \\10^{-8}\\ times the
largest, the repair
[`optimizers7::newton()`](https://statmodels7.github.io/optimizers7/reference/newton.html)
makes, which keeps the increment a descent direction.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
its caller.
