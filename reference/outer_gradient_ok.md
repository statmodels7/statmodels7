# Can the Exact Derivatives of the Criterion Be Computed Here?

`TRUE` when the exact gradient of the marginal criterion in the
hyperparameters under estimation (`order = 1`), or the exact gradient
and Hessian (`order = 2`), can be assembled for the model.

## Usage

``` r
outer_gradient_ok(spec, design, idx, method, order = 1L, gamma = FALSE)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- idx:

  The outer index, from
  [`outer_hyper_index()`](https://statmodels7.github.io/statmodels7/reference/outer_hyper_index.md).

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- order:

  `1` for the gradient, `2` for the Hessian as well.

- gamma:

  `TRUE` where the outer vector also carries coefficients that the
  criterion estimates
  ([`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)),
  so an empty index is not a reason to answer `FALSE`.

## Value

`TRUE` or `FALSE`.

## Details

The family must carry the derivatives that the order needs
([`order_available()`](https://statmodels7.github.io/statmodels7/reference/order_available.md)):
the third derivative of the log-likelihood for the gradient and the
fourth for the Hessian.

**The information.** \\K\\ enters the criterion through its determinant,
so the gradient needs \\\partial K/\partial\beta\\. With the observed
information that is the third derivative of the log-likelihood in the
link-scale predictors. With the expected information it is the
derivative in the predictors of \\-E\[\ell''\]\\, which is not
\\-E\[\ell'''\]\\. That derivative is supplied by
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html),
and the expected route is admitted when the family returns it
([`expected_deriv_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv_ok.md)).
At order 2 the expected route is admitted only for
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
for a model without a block that moves with its coefficients and without
a structural term, and for a family that returns
[`distributions7::distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.html)
([`expected_deriv2_ok()`](https://statmodels7.github.io/statmodels7/reference/expected_deriv2_ok.md)).
The expected route is not admitted at either order where a penalty
covers the parameters of a structural term, because the criterion is
then assembled from the observed joint curvature.

**The structural terms.** A term of the likelihood shape (`regime()` and
the marginal break-point term) has no exact derivatives at either order.
For a term of the filter shape, order 2 requires a penalty over the
term's own parameters and a term that implements
[`modelterms7::term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.html).
Where a penalty covers the term's parameters, the term must also
implement
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html).

**The moving blocks.** A block that moves with its coefficients (`nl()`,
`seg()`, `jump()` and `jseg()`) is admitted at order 1 on either route
and at order 2 on the observed route. Its contribution is read from the
term
([`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md),
[`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md))
and is never differenced.

**The penalties.** \\\partial S/\partial\theta\\ and its second
derivative are given by the generics of penalties7
([`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html),
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html)
and
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)).
The function calls the derivative methods of a penalty at probe values
and reads an error as a penalty that does not supply them, so that a
penalty is admitted according to what its methods return and not
according to its class. The check has two parts. Every penalized unit of
the model is examined, because the movement of the mode reaches the
determinant through each penalty whose Hessian moves with the
coefficients (a Student t prior, for instance). A penalty that
[`penalties7::beta_quadratic()`](https://statmodels7.github.io/penalties7/reference/beta_quadratic.html)
calls quadratic passes, and so does a penalty with a kink (a lasso, a
SCAD or an MCP). Any other penalty must supply
[`penalties7::penalty_dhessian_beta()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian_beta.html),
and a mixed covariance class of this kind is not admitted. Beside a
structural term such a penalty is not admitted at order 2. The second
part examines, with
[`penalty_answers()`](https://statmodels7.github.io/statmodels7/reference/penalty_answers.md),
the penalty of every term that owns a hyperparameter under estimation,
and for a hyperparameter shared by several terms the penalty of each of
them. A penalty with a kink does not remove the exact derivatives of the
other hyperparameters, because its own hyperparameters are left out of
the outer index.

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
