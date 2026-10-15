# One Group's Gaussian Prior From Several Penalties

Assembles the precision of one group's effect where a random-effect term
declares several penalties, each over some of the within-group columns:
one per column under `correlated = FALSE`, or one over the null space
and one over the penalized coordinates of a smooth random effect.

## Usage

``` r
random_prior_joint(units, fit, param, key, d)
```

## Arguments

- units:

  The term's penalized units, from
  [`statmod_penalized()`](https://statmodels7.github.io/statmodels7/reference/statmod_penalized.md).

- fit:

  The fitted model, for the hyperparameters.

- param, key:

  The term.

- d:

  The number of within-group columns.

## Value

A list with `gaussian = TRUE`, `dim` and `chol`, as
[`random_prior()`](https://statmodels7.github.io/statmodels7/reference/random_prior.md)
returns for a Gaussian prior.

## Details

Every penalty of such a term reads the same within-group columns in
every group, group by group, so the first group's positions of each
entry, and the corresponding block of the penalty's Hessian, are that
entry's part of one group's precision \\\Omega\\. The parts cover
disjoint columns and are added. A term whose penalties are not all
Gaussian has no closed-form joint prior and signals an error.
