# How a Moving Block Enters the Twice-Contracted Fourth Derivative

The part of \\\mathrm{tr}(M\\U\[v,u\])\\ that
[`trace_design_form`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
does not compute, where \\U\\ is the second derivative of \\K\\ in the
coefficients contracted in two directions.

## Usage

``` r
trace_refresh4(spec, M, params, Hl, dv, du, units)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- M:

  The matrix the trace is taken against.

- params:

  The parameter names.

- Hl:

  The link-scale curvature, from
  [`refresh_hessian`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md).

- dv, du:

  The two directions, from
  [`refresh_direction`](https://statmodels7.github.io/statmodels7/reference/refresh_direction.md).

- units:

  The refreshable terms, from
  [`refresh_units`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

## Value

A single number.

## Details

Differentiating the three contributions of \\\partial K/\partial\beta\\
once more leaves three further terms that read \\\partial
X/\partial\beta\\ and one that would read
\\\partial^2X/\partial\beta^2\\. Writing \\D_a = (\partial
X_a/\partial\beta)v\\ and \\E_a = (\partial X_a/\partial\beta)u\\, \$\$N
= N_1 + N_2 + N_3, \qquad N_1\[(a,j),(b,k)\] = -\sum_i w_i\Big(\sum_m
\ell\_{abm}(X_mv_m)\_i\Big) E_a\[i,j\]X_b\[i,k\],\$\$ \\N_2\\ the same
with the two directions exchanged, and \\N_3\[(a,j),(b,k)\] = -\sum_i
w_i\ell\_{ab,i}D_a\[i,j\]E_b\[i,k\]\\, which is supported where BOTH
blocks move. The correction is \\\mathrm{tr}(M(N + N^\top)) = 2\sum
M\odot N\\, and the first two pieces are read off
[`refresh_amat`](https://statmodels7.github.io/statmodels7/reference/refresh_amat.md)
with the third derivative in place of the curvature, so neither \\N\\
nor any contraction is assembled.

**The omitted term** carries \\\partial^2X/\partial\beta^2\\, the term's
own third derivative, which no contract supplies. Measured on a bilinear
\\f\\, where it is exactly zero, the corrected Hessian reaches the
reference's own floor.

## See also

[`trace_design_form`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md),
[`contract3_refresh`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md)
