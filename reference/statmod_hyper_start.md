# The Hyperparameters a Specification Starts From

Builds the hyperparameter structure a fit begins from: one entry per
penalized term, each at the midpoint of its own bounds. That is the
probe rule modelterms7 already uses when it reads a penalty's kinks, so
the two layers start a penalty at the same place.

## Usage

``` r
statmod_hyper_start(spec, design = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md),
  whose terms are walked equation by equation.

## Value

A named list with one entry per distribution parameter, each a named
list of numeric vectors keyed by term key, each vector named by that
penalty's own hyperparameters. An equation with no penalized term holds
an empty list.

## Details

A hyperparameter the term holds keeps the held value instead of the
probe. Which ones those are is said by the term and by nothing else,
through
[`modelterms7::term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.html).

The probe is a placeholder. Every hyperparameter left `NULL` on its term
is estimated afterwards, by a marginal criterion or along a path.

## See also

[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
for one penalty's midpoints,
[`statmod_hyper_merge()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_merge.md)
for a caller's overrides,
[`statmod_held()`](https://statmodels7.github.io/statmodels7/reference/statmod_held.md)
for which are held.
