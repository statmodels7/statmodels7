# How a Moving Block Enters the Twice-Contracted Fourth Derivative

Returns the part of \\\mathrm{tr}(M\\U\[v,u\])\\ that
[`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
does not compute, where \\U\\ is the second derivative of \\K\\ in the
coefficients contracted in two directions.

## Usage

``` r
trace_refresh4(spec, M, params, npar, Hl, dv, du, units, G, d3, v, f2, acurv)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- M:

  The matrix against which the trace is taken.

- params, npar:

  The parameter names and block sizes.

- Hl:

  The link-scale curvature, from
  [`refresh_hessian()`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md).

- dv, du:

  The two directions, from
  [`refresh_direction()`](https://statmodels7.github.io/statmodels7/reference/refresh_direction.md).

- units:

  The refreshable terms, from
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- G:

  The leverage diagonal, from
  [`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md).

- d3:

  The third derivatives on the link scale.

- v:

  The first direction over the stacked coefficients.

- f2:

  The block's second derivative in the two directions, one matrix per
  unit, from
  [`refresh_dblock2()`](https://statmodels7.github.io/statmodels7/reference/refresh_dblock2.md).

- acurv:

  The curvature weights, from
  [`refresh_curv_amat()`](https://statmodels7.github.io/statmodels7/reference/refresh_curv_amat.md).

## Value

A single number.

## Details

Differentiating the three contributions of \\\partial K/\partial\beta\\
once more leaves five terms. With \\D_a = (\partial
X_a/\partial\beta)v\\ and \\E_a = (\partial X_a/\partial\beta)u\\, three
of them read only the first derivative, \$\$N_1\[(a,j),(b,k)\] = -\sum_i
w_i\Big(\sum_m \ell\_{abm}(X_mv_m)\_i\Big) E_a\[i,j\]X_b\[i,k\],\$\$
\\N_2\\ is the same expression when the two directions are exchanged,
and \\N_3\[(a,j),(b,k)\] = -\sum_i
w_i\ell\_{ab,i}D_a\[i,j\]E_b\[i,k\]\\, which is supported where both
blocks move. A fourth reads the second, \\N_4\[(a,j),(b,k)\] = -\sum_i
w_i\ell\_{ab,i}F_a\[i,j\]X_b\[i,k\]\\ with \\F_a =
(\partial^2X_a/\partial\beta^2)\[v,u\]\\ from
[`modelterms7::term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.html).
Those four enter as \\\mathrm{tr}(M(N + N^\top)) = 2\sum M\odot N\\.
\\N_1\\, \\N_2\\ and \\N_4\\ are read off
[`refresh_amat()`](https://statmodels7.github.io/statmodels7/reference/refresh_amat.md),
with the third derivative for the first two and the curvature for the
fourth, and \\N_3\\ is summed from the two blocks' derivatives with the
rows of \\M\\ between them, so that the function assembles neither \\N\\
nor any contraction.

The fifth is of another shape and is added separately. The
third-derivative contraction carries the direction's predictor
\\(X_mv_m)\_i\\, which moves with \\\beta\\ as everything else does
where \\X_m\\ is a Jacobian, and differentiating it leaves \$\$-\sum_i
w_i\Big(\sum_m \ell\_{abm}\\\dot v_m(i)\Big)X_a\[i,j\]X_b\[i,k\], \qquad
\dot v_m = (\partial X_m/\partial\beta\[u\])\\v_m.\$\$ That is the shape
that
[`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
already computes, so the piece is one further call of it with \\\dot v\\
in place of the predictors. \\\dot v\\ is the predictor's own second
derivative contracted in both directions and is therefore symmetric in
them, which is a free check on the assembly.

## See also

[`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md),
[`contract3_refresh()`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md),
[`refresh_dblock2()`](https://statmodels7.github.io/statmodels7/reference/refresh_dblock2.md)
