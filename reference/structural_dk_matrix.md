# The Joint Derivative of the Penalized Information Along One Direction

\\\partial K/\partial u\\\[v\]\\ assembled as a matrix over the
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

The gradient needs the trace and nothing else, so it never forms this.
The HESSIAN needs the matrix twice over: in \\\mathrm{tr}(M K_l M
K_m)\\, which no contraction reduces, and as the operator carrying the
mode's second movement. Writing \\V_a\\ for each equation's rows,
\\\mathrm{d}\varphi = E v\\ for the filter's own moving, and \\E\\ for
the second derivative of the predictor the filter produces,

\$\$\frac{\partial K}{\partial u}\[v\] = -\sum_i
w_i\Big\[\sum\_{a,b}\Big(\sum_k\ell\_{abk}(V_k\cdot v)\Big) V_a^\top
V_b + \sum_b \ell\_{pb}\big(\mathrm{d}\varphi\otimes V_b +
V_b\otimes\mathrm{d}\varphi\big)\Big\] - W(\kappa_v) - W_3\[v\],\$\$

with \\\kappa_v = \sum_k \ell\_{pk}(V_k\cdot v)\\ the weight the level's
own term is re-read at and \\W_3\[v\]\\ what
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
returns. Traced against \\M\\ it reproduces
[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md)
exactly, which is what a test asserts of it.

## See also

[`structural_chain_extra()`](https://statmodels7.github.io/statmodels7/reference/structural_chain_extra.md),
[`statmod_structural_hess()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_hess.md)
