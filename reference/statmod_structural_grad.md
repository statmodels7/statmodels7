# The Exact Gradient Where a Penalty Covers a Filter's Own Parameters

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md)
over the joint vector of coefficients and a structural term's
parameters, over which the determinant spans there.

## Usage

``` r
statmod_structural_grad(
  spec,
  design,
  coef,
  hyper,
  method,
  idx,
  basis = NULL,
  free = TRUE
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients at the penalized mode.

- hyper:

  The hyperparameters.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- idx:

  The outer index.

- basis:

  The integrated subspace, or `NULL`.

- free:

  Whether to carry the result onto the free scale.

## Value

A numeric vector, one entry per row of `idx`, or `NULL` where the
determinant does not exist.

## Details

The criterion is the same and so are its three pieces; what changes is
that \\K\\ carries the recursion's own second derivative, \$\$K =
-\sum_t w_t\sum\_{a,b}\ell\_{ab,t}V\_{a,t}^\top V\_{b,t} - \sum_t w_t
\ell\_{p,t}E_t,\$\$ so differentiating it along the direction \\v\\ the
mode moves in gives three contributions instead of one:

1.  the family's third derivative against the per-observation diagonal
    of \\M\\, which is
    [`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md)'s
    formula with \\V\\ in place of \\X\\;

2.  the derivative of \\V_p\\ itself, which is \\E_t v\\, one row per
    observation, entering through every \\\ell\_{pb}\\;

3.  the derivative of \\\ell_p E_t\\, whose first half is \\E\\
    re-weighted by \\\sum_k \ell\_{pk}(V_k\cdot v)\\ and whose second is
    \\\partial^3 e_t/\partial u^3\[v\]\\, the one genuinely new object.

Piece 1 is direction-free and is computed once; pieces 2 and 3 are read
along each hyperparameter's own direction, which is the trade that keeps
the cost at \\O(nm^2)\\ per hyperparameter instead of \\O(nm^3)\\ once.

## See also

[`statmod_marginal_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_grad.md),
[`modelterms7::term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.html)
