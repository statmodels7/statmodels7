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
  A simulation drawing a prior over coordinates that ride a chart passes
  a half: measured, a scale of 1 puts 10.1 per cent of the persistences
  a partial-autocorrelation chart carries past 0.95 and a scale of 0.5
  puts 0.9 per cent there.

## Value

A named numeric vector, one entry per hyperparameter of `pen`, named as
the penalty names them. `numeric(0)` for a penalty with none, which a
fixed prior is.

## Details

One is the scale a smoothing parameter lives on before anything is known
about it, and that is an argument about a hyperparameter bounded below,
which lives on \\(0, \infty)\\ and for which one is a scale. A
coordinate unbounded on both sides is not such a thing: it is a free
coordinate of a parameters7 chart, and there one is not one. On
`dr_prod(2)` it reads as standard deviations of \\e^1 = 2.718\\ and a
correlation of -0.664, where zero is the chart's neutral point – unit
standard deviations and no correlation – which is where an ordinary
random effect starts, its own scale being bounded below and starting at
`lower + 1`.

The value matters only as somewhere to begin: it is a probe, and the
criterion moves it at the first opportunity. What it can do is leave the
inner fit somewhere it cannot converge from, and measured it did: over
22 models across the four charts a caller can reach – `log_cholesky`,
`dr_prod`, `ar1` and `compound_symmetry` – the two starting points give
criteria agreeing to 3e-05 or better wherever both finish, and three
models finish from the neutral point that do not finish from one. A
hyperparameter with a finite bound is untouched, so a ridge, a smooth
and an ordinary random effect all begin exactly where they did.

## See also

[`statmod_hyper_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_hyper_start.md),
its caller.
