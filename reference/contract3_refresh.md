# How a Moving Block Enters the Contracted Third Derivative

The part of \\T\[v\] = (\partial K/\partial\beta)\cdot v\\ that
[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
does not compute, as a matrix over the stacked coefficients.

## Usage

``` r
contract3_refresh(spec, design, params, npar, offs, total, dir, Hl, units)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- params, npar, offs, total:

  The block bookkeeping.

- dir:

  One entry per unit, from
  [`refresh_direction()`](https://statmodels7.github.io/statmodels7/reference/refresh_direction.md).

- Hl:

  The link-scale curvature, from
  [`refresh_hessian()`](https://statmodels7.github.io/statmodels7/reference/refresh_hessian.md).

- units:

  The refreshable terms, from
  [`refresh_units()`](https://statmodels7.github.io/statmodels7/reference/refresh_units.md).

## Value

A square matrix, whose transpose completes the correction.

## Details

Differentiating \\H\\ in \\\beta\\ and contracting with \\v\\ leaves,
beside the third derivative
[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md)
carries, two terms in \\\partial X/\partial\beta\\: with \\D_a =
(\partial X_a/\partial\beta)v\\, \$\$R\[(a,j),(b,k)\] = -\sum_i
w_i\\\ell\_{ab,i}\\D_a\[i,j\]\\X_b\[i,k\],\$\$ and the other is its
transpose, so the correction is \\R + R^\top\\. The derivative is asked
of the term through
[`modelterms7::term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.html)
and never differenced here, for the reason
[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
records: a break-point column is a step function in its break-point.

## See also

[`contract3()`](https://statmodels7.github.io/statmodels7/reference/contract3.md),
[`u_refresh()`](https://statmodels7.github.io/statmodels7/reference/u_refresh.md)
