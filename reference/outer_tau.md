# The Effective Degrees of Freedom of a Whole Fit

Computes \\\tau = \mathrm{tr}\[(H+S)^{-1}H\]\\, the trace of the
influence matrix, which is the number of effective degrees of freedom
charged by a prediction-error criterion. An unpenalized fit has \\S =
0\\ and \\\tau = p\\; shrinkage lowers it, and \\\tau\\ falls to the
dimension of the penalty's null space as the hyperparameter grows.

## Usage

``` r
outer_tau(J, H, active = NULL)
```

## Arguments

- J:

  The penalized information \\H + S\\, a `p x p` symmetric matrix, dense
  or sparse.

- H:

  The likelihood's own information, of the same shape as `J`.

- active:

  A logical vector over the stacked coefficients marking the ones away
  from a kink, or `NULL` (the default) for all of them. Both matrices
  are subset to it before the trace.

## Value

A single number: the trace, `0` when `active` selects nothing, and
`NA_real_` when `J` has no Cholesky factor.

## Restricting to the active set

A penalty with a kink is not twice differentiable everywhere, and the
mode of such a fit sits on the kink for every coefficient that it sets
to zero. It is twice differentiable away from the kink, so the trace is
taken over the active coordinates only, which form the submodel that the
fit has selected.

For a penalty that is linear there, the lasso, \\S\_{AA}\\ vanishes and
\\\tau\\ is exactly the number of coefficients that survived, the
unbiased count of Zou, Hastie and Tibshirani (2007). The elastic net
keeps its quadratic part and spends slightly less; SCAD and MCP
contribute their own curvature and spend slightly more.

The trace over the whole vector does not depend on which coefficients
are at the kink, so a criterion built on it assigns the same degrees of
freedom to a model that selects nothing and to one that selects
everything.

## The inverse

The inverse is computed through a Cholesky factorization of \\J\\, and
the result is `NA_real_` when that fails, which a caller reads as an
unusable point. Both matrices are densified first: the trace needs the
full inverse, which is dense whatever \\J\\ was.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom
of the lasso. *The Annals of Statistics* 35(5), 2173–2192.

## See also

[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
for the criterion that this function feeds,
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
for the per-term reading, which splits the same trace by term.
