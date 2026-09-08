# A Covariance Class Split Between the Coefficients and a Filter's Parameters

The value, gradient and Hessian of every penalty whose coordinates are
part coefficients and part a structural term's own parameters, placed in
the joint vector \\\[\beta; \zeta\_{\mathrm{free}}\]\\.

## Usage

``` r
joint_penalty_at(
  spec,
  design,
  coef,
  hyper,
  what = c("value", "gradient", "hessian"),
  n = 0L
)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  A named list of coefficient vectors.

- hyper:

  The hyperparameters.

- what:

  Which quantity: `"value"`, `"gradient"` or `"hessian"`.

- n:

  The length of the joint vector, `nb` plus the number of free
  parameters. Read from the caller, which has already computed it.

## Value

For `"value"` a single number; for `"gradient"` a numeric vector of
length `n`; for `"hessian"` an `n` by `n` matrix. Zero throughout where
no class is mixed.

## Why one function and not two halves

A class's prior is one object over one stacked vector. penalties7
computes one gradient and one Hessian over that vector, in the class's
own order, and knows nothing of the split – to it these are \\m\\ blocks
of \\d\\ numbers with a covariance in common. Splitting the answer into
a coefficient half and a parameter half would report each as though the
other were not there, and would lose the CROSS block entirely, which is
the only place the correlation between the two enters at all.

So the scatter is written once and the three assemblies that build the
joint matrix add what it returns: the inner step
([`statmod_fit_joint()`](https://statmodels7.github.io/statmodels7/reference/statmod_fit_joint.md)),
the marginal criterion
([`statmod_marginal_full()`](https://statmodels7.github.io/statmodels7/reference/statmod_marginal_full.md))
and the variance
([`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)).
[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
and
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
skip a mixed class for the same reason.

## Where the values come from

Each coordinate is read where it lives – a coefficient from `coef`, a
filter parameter from the design's structural state – and assembled in
the class's own order, which `joint` on the unit already records. The
result is written back at those same positions, so no ordering is
composed twice.

## See also

[`statmod_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalty_at.md)
for the coefficients' own penalties,
[`statmod_structural_penalty()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_penalty.md)
for a filter's.
