# The Joint Derivative of the Penalized Information Along One Direction

Returns \\\partial K/\partial u\\\[v\]\\ assembled as a matrix over the
coefficients and the filter's own parameters, where
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
returns only its trace against \\M\\.

## Usage

``` r
structural_dk_matrix(spec, design, jd, st, v)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- jd:

  The joint rows, from
  [`joint_design_rows()`](https://statmodels7.github.io/statmodels7/reference/joint_design_rows.md).

- st:

  The shared quantities, from
  [`structural_grad_parts()`](https://statmodels7.github.io/statmodels7/reference/structural_grad_parts.md).

- v:

  The direction, over the estimated coordinates.

## Value

A square symmetric matrix over the estimated coordinates.

## Details

The gradient uses only the trace of this matrix against \\M\\ and does
not form the matrix. The Hessian uses the matrix in two places: in
\\\mathrm{tr}(M K_l M K_m)\\, which a contraction cannot reduce, and as
the operator that carries the second movement of the mode. Write \\V_a\\
for the rows of each equation, \\\mathrm{d}\varphi = E v\\ for the
movement of the filter's own parameters, and \\E\\ for the second
derivative of the predictor that the filter produces. Then

\$\$\frac{\partial K}{\partial u}\[v\] = -\sum_i
w_i\Big\[\sum\_{a,b}\Big(\sum_k\ell\_{abk}(V_k\cdot v)\Big) V_a^\top
V_b + \sum_b \ell\_{pb}\big(\mathrm{d}\varphi\otimes V_b +
V_b\otimes\mathrm{d}\varphi\big)\Big\] - W(\kappa_v) - W_3\[v\],\$\$

with \\\kappa_v = \sum_k \ell\_{pk}(V_k\cdot v)\\ the weight at which
the term of the level is read again and \\W_3\[v\]\\ the curvature that
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
returns. The trace of the matrix against \\M\\ equals the value of
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
up to rounding.

## See also

[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md),
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
