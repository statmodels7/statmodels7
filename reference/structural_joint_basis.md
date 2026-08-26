# The Subspace a Marginal Criterion Integrates Over, Jointly

The basis
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
projects the joint curvature onto: the coefficients' own range basis,
with one column per penalized parameter of the structural term appended.

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

A matrix.

## Details

It is written once because
[`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md)
and
[`statmod_structural_grad()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_grad.md)
must project onto the same subspace, and two callers composing it
separately would agree only by accident.
