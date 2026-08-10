# Can the Exact Gradient Be Computed Here?

`TRUE` when every hyperparameter under estimation belongs to a penalty
whose second derivative in the coefficients is linear in the
hyperparameters and free of the coefficients, and the criterion uses the
observed information.

## Usage

``` r
outer_gradient_ok(spec, design, idx, method)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- idx:

  The outer index, from
  [`outer_hyper_index`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

## Value

`TRUE` or `FALSE`.

## Details

**Why the observed information.** \\K\\ enters the criterion through its
determinant, so the gradient needs \\\partial K/\partial\beta\\. With
the observed information that is the third derivative of the
log-likelihood in the link-scale predictors, which every family of
distributions7 carries in closed form. With the expected information it
would be the derivative in \\\beta\\ of \\-E\[\ell''\]\\, which is not
\\-E\[\ell'''\]\\ and is not one of that package's generics. So
`reml(hessian = "observed")` is what the exact route asks for, and
`"expected"` keeps the derivative-free search.

**Why linear.** \\\partial S/\partial\theta\\ is not a generic of
penalties7: what that package exposes is the penalty, its gradient, its
Hessian and the mixed block, not the third derivative
\\\partial^3\rho/\partial\beta^2\partial\theta\\. For a penalty whose
Hessian is linear in the hyperparameters the derivative is recoverable
from the Hessian itself and nothing has to be differentiated; for one
that is not – a ridge, whose Hessian carries \\1/\sigma^2\\, or any
penalty built from a density – it is not, and the search stays
derivative-free. Closing that case is a generic in penalties7, and for
the density branch it would need \\\partial^3\ell/\partial
y^2\partial\theta\\ from distributions7, which does not exist either.

The linearity is **checked and not assumed**: a penalty is asked for its
Hessian at \\\theta\\ and at \\2\theta\\, and at two coefficient
vectors, and admitted only if the first doubles and the second does not
move.

## See also

[`statmod_marginal_grad`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
