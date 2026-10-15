# A Penalty's Starting Hyperparameters

Picks a starting value for each of a penalty's hyperparameters from its
`params_bounds`, one unit inside whichever ends are finite:

- both ends finite: their midpoint;

- bounded below only, the common case: `lower + unit`, so a
  hyperparameter on \\\[0, \infty)\\ starts at 1 at the default `unit`;

- bounded above only: `upper - unit`;

- unbounded on both sides: 0, whatever `unit` is.

## Usage

``` r
penalty_theta_start(pen, unit = 1)
```

## Arguments

- pen:

  A penalties7 penalty object, read for its `params_bounds` property
  alone.

- unit:

  How far from a one-sided bound the neutral point sits, `1` by default.
  A simulation that draws a prior over the coordinates of a chart may
  pass a smaller value, which keeps the persistences drawn on a
  partial-autocorrelation chart away from the edge of the chart.

## Value

A named numeric vector, one entry per hyperparameter of `pen`, named as
the penalty names them. `numeric(0)` for a penalty with none, which a
fixed prior is.

## Details

A hyperparameter bounded below lives on \\(0, \infty)\\ when the lower
bound is zero, and a start of one is a scale on that range. A coordinate
that is unbounded on both sides is a free coordinate of a parameters7
chart, and for a chart the value one has no such meaning: on
`dr_prod(2)` it gives standard deviations of \\e^1 = 2.718\\ and a
correlation of -0.664. Zero is the neutral point of the chart, with unit
standard deviations and no correlation, so that a coordinate unbounded
on both sides starts at zero. This is also where an ordinary random
effect starts, because its own scale is bounded below and starts at
`lower + 1`.

The value is a starting point only, and the criterion moves it. A start
that the inner fit cannot converge from can end the search, and the
neutral point of an unbounded coordinate converges in cases where a
start at one does not. A hyperparameter bounded on both sides starts at
the midpoint of its bounds, and one bounded on one side starts one
`unit` inside the bound. The smoothing parameter of a smooth, the
precision of a ridge, the penalty of a lasso and the scale of an
ordinary random effect are bounded below by zero and start at one.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md),
its caller.
