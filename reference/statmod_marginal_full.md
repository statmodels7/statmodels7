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
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients, at the penalized mode.

- hyper:

  The hyperparameters.

- basis:

  The integrated subspace over the coefficients, or `NULL` for
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md).

## Value

A square matrix, or `NULL` where the term could not be run.

## Details

A marginal criterion integrates the quantities its penalty shrinks.
Where that penalty is a term's own – the deviations of a panel – those
quantities are not coefficients, so the determinant has to span them
too: taken over the coefficients alone it does not depend on the
hyperparameter at all, and the criterion is the penalized likelihood,
whose maximum in a shrinkage parameter is at no shrinkage.

Both pieces exist already.
[`statmod_full_information`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
spans the coefficients followed by the term's free parameters, which is
the order the joint fit uses, and
[`statmod_structural_penalty`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
gives the penalty's Hessian in those same parameters. The information
here is the OBSERVED one whatever the method asks for: the expected
information over a filter's own parameters is not one of the quantities
the toolkit carries, the recursion's state entering the expectation.

For
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
the tail is integrated over the penalized coordinates alone, through the
penalty's own range basis, so a parameter the penalty does not cover – a
population value – is profiled rather than integrated, exactly as an
unpenalized coefficient is.

## See also

[`statmod_marginal`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal.md)
