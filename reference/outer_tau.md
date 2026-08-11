# The Effective Degrees of Freedom of a Whole Fit

\\\tau = \mathrm{tr}\[(H+S)^{-1}H\]\\, the trace of the influence
matrix.

A penalty with a kink is not twice differentiable everywhere, and the
mode of such a fit sits at the kink for every coefficient it sets to
zero. It is twice differentiable away from the kink, so the trace is
taken over the active coordinates alone, which is the submodel the fit
has selected. For a penalty that is linear there – the lasso –
\\S\_{AA}\\ vanishes and \\\tau\\ is the number of coefficients that
survived, which is the unbiased count of Zou, Hastie and Tibshirani
(2007). The elastic net keeps its quadratic part, and SCAD and MCP their
own curvature.

## Usage

``` r
outer_tau(J, H, active = NULL)
```

## Arguments

- J:

  The penalized information.

- H:

  The likelihood's information.

- active:

  A logical vector over the stacked coefficients, or `NULL` for all of
  them.

## Value

A single number.

## References

Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of freedom
of the lasso. *The Annals of Statistics* 35(5), 2173–2192.
