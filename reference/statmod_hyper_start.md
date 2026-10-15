# The Hyperparameters a Specification Starts From

Builds the hyperparameter structure that a fit starts from: one entry
per penalized term, each hyperparameter at the value given by
[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md).
That value is the probe value that modelterms7 uses when it reads the
kinks of a penalty, so that both layers start a penalty at the same
place.

## Usage

``` r
statmod_hyper_start(spec, design = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  whose terms are walked equation by equation.

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it. It is built from `spec` when `NULL`.

## Value

A named list with one entry per distribution parameter, each a named
list of numeric vectors keyed by term key, each vector named by that
penalty's own hyperparameters. An equation with no penalized term holds
an empty list.

## Details

A hyperparameter that its term holds takes the value that the term
carries. A value held on one member of a group of terms that share a
hyperparameter through a common label applies to every member of the
group.

The start value is a placeholder. Every hyperparameter that is `NULL` on
its term is estimated afterwards, by a marginal criterion or along a
path.

## See also

[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
for the start values of one penalty,
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
for the hyperparameters that are held.
