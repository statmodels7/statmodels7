# The Level Chart of a Filter's Joint Fit

The change of coordinates under which
[`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)
searches a free level: the level \\\omega\\ of a score-driven term is
replaced by \\L = \omega / \prod_k (1 - r_k)\\, the level around which
the predictor fluctuates, \\r_k\\ being the partial autocorrelations
that carry the persistence (by the Durbin-Levinson recursion, \\1 -
\sum_k \phi_k = \prod_k (1 - r_k)\\). The other coordinates are
unchanged.

## Usage

``` r
joint_level_chart(spec, jp)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- jp:

  The joint pieces, as
  [`statmod_joint_pieces()`](https://statmodels7.github.io/statmodels7/reference/statmod_joint_pieces.md)
  returns them.

## Value

`NULL`, or a list of the functions `to` (chart to the term's
coordinates), `from`, `grad` and `hess`.

## Details

The chart applies where the level is one free coordinate on the identity
link and every persistence is a free scalar partial autocorrelation; in
any other case the result is `NULL` and the fit runs on the term's own
coordinates. With \\\omega = L P(\rho)\\, the gradient is \\J^\top g\\
and the Hessian \\J^\top H J + g\_\omega \nabla^2 \omega\\, \\J\\ the
Jacobian of the map.
