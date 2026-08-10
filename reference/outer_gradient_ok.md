# Can the Exact Gradient Be Computed Here?

`TRUE` when every hyperparameter under estimation belongs to a penalty
whose second derivative in the coefficients is linear in the
hyperparameters and free of the coefficients, and the criterion uses the
observed information.

## Usage

``` r
outer_gradient_ok(spec, design, idx, method, order = 1L)
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

- order:

  `1` for the gradient, `2` for the Hessian as well.

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

**Why the penalty is asked.** \\\partial S/\partial\theta\\ and its
second derivative are generics of penalties7
([`penalty_dhessian`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html),
[`penalty_d2hessian`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html),
[`penalty_dcross`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)).
A penalty that answers them is estimable by a marginal criterion
whatever its shape: the quadratic, the additive, the structured and the
separable branches all do, so a ridge, a random effect and a
heavy-tailed prior are covered as well as a spline. One that does not –
a SCAD, an MCP, anything with a kink – rejects, and the search stays
derivative-free. Nothing here tests a penalty's behaviour to find out
what it is; it is asked.

## See also

[`statmod_marginal_grad`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
