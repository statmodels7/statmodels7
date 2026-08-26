# The Penalized Curvature Over the Coefficients and a Filter's Parameters

The matrix whose determinant the Laplace approximation reads, where a
penalty covers a structural term's own parameters instead of a block of
design columns.

## Usage

``` r
statmod_marginal_full(spec, design, coef, hyper, basis = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, at the penalized mode.

- hyper:

  The hyperparameters.

- basis:

  The integrated subspace over the coefficients, as
  [`integrated_basis()`](https://statmodels7.github.io/statmodels7/reference/integrated_basis.md)
  returns it, or `NULL` for
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

## Value

A square symmetric matrix spanning the integrated coefficients followed
by the integrated structural parameters. `NULL` where the structural
term could not be run at the given point, which a caller reads as an
unusable hyperparameter.

## Why the determinant has to span them

A marginal criterion integrates the quantities its penalty shrinks.
Where the penalty is a structural term's own, as with the deviations of
a panel, those quantities are not coefficients.

Taken over the coefficients alone the determinant then does not depend
on the hyperparameter at all, and the criterion reduces to the penalized
likelihood, whose maximum in a shrinkage parameter is at no shrinkage.
It does not merely lose accuracy: it answers a different question.

## Both pieces already exist

[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
spans the coefficients followed by the term's free parameters, which is
the order the joint fit uses, and
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
gives the penalty's Hessian in those same parameters.

The information here is the **observed** one whatever the method asks
for. The expected information over a filter's own parameters is not a
quantity the toolkit carries, the recursion's state entering the
expectation.

## Under ML

The tail is integrated over the penalized coordinates alone, through the
penalty's own range basis, so a parameter the penalty does not cover,
such as a population value, is profiled instead of integrated. That is
what happens to an unpenalized coefficient too.

## See also

[`statmod_marginal()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md),
the caller,
[`structural_penalized()`](https://statmodels7.github.io/statmodels7/reference/structural_penalized.md)
for the predicate that selects this route.
