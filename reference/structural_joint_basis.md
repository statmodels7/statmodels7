# The Joint Subspace over Which a Marginal Criterion Integrates

Returns the basis onto which
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
projects the joint curvature: the coefficients' own range basis, to
which one column per penalized parameter of the structural term is
appended.

## Usage

``` r
structural_joint_basis(spec, design, key, free, nb, basis)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- key:

  The structural term's name.

- free:

  Its free parameters, in order.

- nb:

  The number of coefficients.

- basis:

  The coefficients' basis.

## Value

A matrix with `nb` plus `length(free)` rows.

## Details

It is written once because
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
and
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md)
must project onto the same subspace, and a separate composition in each
of the two callers would agree only by accident.
