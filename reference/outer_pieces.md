# The Per-Hyperparameter Pieces of the Outer Derivatives

Returns the penalty's contributions to the criterion's gradient and
Hessian, placed in the stacked coefficient space.

## Usage

``` r
outer_pieces(spec, design, coef, hyper, idx, offs, total, order = 2L)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- hyper:

  The hyperparameters.

- idx:

  The outer index.

- offs, total:

  The block offsets and the total width.

- order:

  `1` for the first-order pieces only, `2` for the second-order ones as
  well.

## Value

A list with `S` (one matrix per hyperparameter) and `c` (one vector),
and at order 2 also `S2` and `c2` (one per pair, keyed), `rho2` (the
hyperparameter Hessian) and `pair` (the key of each pair).

## Details

At order 1 the function calls
[`penalties7::penalty_dhessian()`](https://statmodels7.github.io/penalties7/reference/penalty_dhessian.html)
and
[`penalties7::penalty_cross()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
only. The second-order generics
[`penalties7::penalty_d2hessian()`](https://statmodels7.github.io/penalties7/reference/penalty_d2hessian.html),
[`penalties7::penalty_dcross()`](https://statmodels7.github.io/penalties7/reference/penalty_dcross.html)
and
[`penalties7::penalty_hess_theta()`](https://statmodels7.github.io/penalties7/reference/penalty_grad_theta.html)
are called at order 2, so a penalty that supplies the first two generics
gives an exact gradient without implementing the others.

A shared hyperparameter occupies one row of `idx` and receives the sum
of the contributions of its members. A penalty over the own parameters
of a structural term and a mixed covariance class contribute nothing
here;
[`structural_outer_pieces()`](https://statmodels7.github.io/statmodels7/reference/structural_outer_pieces.md)
places them on the joint vector.
