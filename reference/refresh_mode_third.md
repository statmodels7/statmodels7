# How a Moving Block Enters the Mode's Second Movement

The part of \\\partial^3L/\partial\beta^3\[v,u\]\\ that
[`contract3`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
and
[`contract3_refresh`](https://statmodels7.github.io/statmodels7/reference/contract3_refresh.md)
do not carry, as a vector over the stacked coefficients.

## Usage

``` r
refresh_mode_third(spec, params, npar, units, Hl, gl, tv, du, f2, total)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- params, npar, total:

  The block bookkeeping.

- units:

  The refreshable terms, from
  [`refresh_units`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

- Hl:

  The link-scale curvature, from
  [`refresh_hessian`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md).

- gl:

  The link-scale score.

- tv:

  The predictors of the first direction.

- du:

  The second direction, from
  [`refresh_direction`](https://statmodels7.github.io/statmodels7/reference/refresh_direction.md).

- f2:

  The block's second derivative, one matrix per unit.

## Value

A numeric vector as long as the stacked coefficients.

## Details

\\b\_{ml}\\ solves \\J b\_{ml} = -\[(S_l + T\[b_l\])b_m + S_m b_l +
c\_{ml}\]\\ with \\T\\ the THIRD derivative of the penalized objective
in \\\beta\\. Where a block moves, the objective's second derivative is
not \\K\\ but \\K + D\\, with \\D\\ the term
[`mode_curvature`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md)
builds, so its third derivative is not \\\partial K/\partial\beta\\
either. What is missing is \\\partial D/\partial\beta\\, and
differentiating \\D\_{cd} = -\sum_i w_i\ell_a\\\Xi_a\[i,c,d\]\\ once and
contracting gives two pieces, \$\$-\sum_i w_i\Big(\sum_k
\ell\_{ak}(X_kv_k)\_i\Big)\\ \big\[(\partial
X_a/\partial\beta)u\big\]\[i,c\] \\-\\\sum_i
w_i\\\ell_a(i)\\F_a\[i,c\],\$\$ the first reading the block's first
derivative and the second its SECOND, \\F_a =
(\partial^2X_a/\partial\beta^2)\[v,u\]\\. Both are one `crossprod`
against a per-observation weight, so neither the third-derivative array
of the predictor nor any contraction of it is formed.

**Measured** against the mode refitted at four hyperparameter values and
differenced twice, on `nl(a ~ 0 + ridge(~grp))`: \\b_m\\ is right to
3.8e-08 without this and \\b\_{ml}\\ is wrong by 7.6 to 9.0 per cent at
every step tried, a systematic error and not the reference's noise. With
it the gap falls to 1.8e-05 at the reference's best step, which is
inside the spread between the reference's own consecutive steps.

## See also

[`mode_curvature`](https://statmodels7.github.io/statmodels7/reference/mode_curvature.md),
[`trace_refresh4`](https://statmodels7.github.io/statmodels7/reference/trace_refresh4.md)
