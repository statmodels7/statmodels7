# The Exact Gradient of the Marginal Criterion

\\\partial V/\partial\eta\\ at the penalized mode, over the free scale
of the hyperparameters under estimation.

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
  ctx = NULL
)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The outer index.

- basis:

  The integrated subspace, or `NULL`.

- free:

  Whether to carry the result onto the free scale. The Hessian asks for
  the parameter scale, having its own second-order chain rule to apply.

## Value

A numeric vector, one entry per row of `idx`, or `NULL` where the
determinant does not exist.

## Details

The three pieces are the envelope term \\-\partial\rho/\partial\theta\\,
the explicit derivative of the determinant \\\mathrm{tr}(M\\\partial
S/\partial\theta)\\, and the implicit one \\u'v\\, where \\v =
-(H+S)^{-1}\partial^2\rho/\partial\beta\partial\theta\\ is how the mode
moves and \\u_c = \mathrm{tr}(M\\\partial K/\partial\beta_c)\\ is how
the determinant reads that movement.

\\u\\ is assembled without forming any third-derivative array. Writing
\\G\_{ab,i} = x\_{ia}'M\_{\[a\]\[b\]}x\_{ib}\\ for the per-observation
diagonal of the block of \\M\\, \$\$u\_{k} = -X_k'\Big(w \sum\_{a,b}
\ell'''\_{abk}\\ G\_{ab}\Big),\$\$ one crossprod per distribution
parameter. The component \\\ell'''\_{abk}\\ is looked up by a name BUILT
from the parameter names in the family's own order, never parsed out of
one.

## See also

[`statmod_marginal`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
[`reml`](https://statmodels7.github.io/statmodels7/reference/reml.md)
