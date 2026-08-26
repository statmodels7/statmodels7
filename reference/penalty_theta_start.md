# A Penalty's Starting Hyperparameters

Picks a starting value for each of a penalty's hyperparameters from its
`params_bounds`, one unit inside whichever ends are infinite:

- both ends finite: their midpoint;

- bounded below only, the common case: `lower + 1`, so a hyperparameter
  on \\\[0, \infty)\\ starts at 1;

- bounded above only: `upper - 1`;

- unbounded: 1.

## Usage

``` r
penalty_theta_start(pen)
```

## Arguments

- pen:

  A penalties7 penalty object, read for its `params_bounds` property
  alone.

## Value

A named numeric vector, one entry per hyperparameter of `pen`, named as
the penalty names them. `numeric(0)` for a penalty with none, which a
fixed prior is.

## Details

One is the scale a smoothing parameter lives on before anything is known
about it. The value matters only as somewhere to begin: it is a probe,
and the criterion moves it at the first opportunity.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md),
its caller.
