# The Exact Gradient of the Marginal Criterion

Computes \\\partial V/\partial\eta\\ at the penalized mode, over the
free scale of the hyperparameters under estimation.

## Usage

``` r
statmod_marginal_grad(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  free = TRUE,
  ctx = NULL,
  gam = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The outer index.

- basis:

  The integrated subspace, or `NULL`.

- free:

  Whether to carry the result onto the free scale. The Hessian assembly
  calls this function with `free = FALSE` and applies its own
  second-order chain rule.

- ctx:

  The evaluation context, or `NULL`.

- gam:

  The coefficients that the criterion estimates, from
  [`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md),
  or `NULL`.

## Value

A numeric vector, one entry per row of `idx` followed by one per
position of `gam$where`, or `NULL` where the determinant does not exist.

## Details

The gradient is the sum of the envelope term
\\-\partial\rho/\partial\theta\\, the explicit derivative of the
determinant, \\-\tfrac12\mathrm{tr}(M\\\partial S/\partial\theta)\\, the
implicit one, \\-\tfrac12 u'v\\, and, where a penalty's Hessian depends
on the coefficients, \\-\tfrac12\mathrm{tr}(M\\\partial
S/\partial\beta\[v\])\\. Here \\v =
-J^{-1}\\\partial^2\rho/\partial\beta\partial\theta\\ is how the mode
moves, \\u_c = \mathrm{tr}(M\\\partial K/\partial\beta_c)\\ is how the
determinant reads that movement, and \\M\\ is \\K^{-1}\\ under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and the inverse of \\K\\ on the integrated subspace under
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

The matrix \\J\\ is the curvature of the penalized likelihood at the
mode: the observed information plus the penalty's Hessian, plus the
second derivative of a block that moves with its coefficients
([`mode_curvature()`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md)).
It differs from \\K\\ when the criterion reads the expected information.
A coefficient that the criterion estimates, a direction that the model
does not identify and a coordinate that a kink of a break-point term
holds do not move with \\\theta\\, and the function solves on the
directions that they leave free
([`mode_solver()`](https://statmodels7.github.io/statmodels7/reference/mode_solver.md)).
Where a penalty covers the parameters of a structural term, the function
delegates to
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md).

\\u\\ is assembled without forming any third-derivative array. Writing
\\G\_{ab,i} = x\_{ia}'M\_{\[a\]\[b\]}x\_{ib}\\ for the per-observation
diagonal of the block of \\M\\, \$\$u\_{k} = -X_k'\Big(w \sum\_{a,b}
\ell'''\_{abk}\\ G\_{ab}\Big),\$\$ one crossprod per distribution
parameter. The component \\\ell'''\_{abk}\\ is looked up by a name built
from the parameter names in the family's own order, never parsed out of
one.

For each coefficient in `gam`, which the criterion estimates, the entry
is the score of the log-likelihood in that coefficient minus half of the
movement of the determinant along the direction in which the mode moves
when the coefficient changes by one unit.

## See also

[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
