# Draw One Penalty's Hyperparameters

A value per hyperparameter, drawn on the chart that the penalty carries
for it around the neutral point that
[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
gives.

## Usage

``` r
rstatmod_hyper(pen, sd, unit = 1)
```

## Arguments

- pen:

  A penalties7 penalty.

- sd:

  The width of the draws.

- unit:

  How far from a one-sided bound the draw is centered, passed to
  [`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md).
  It is one half for a prior over a structural term's own parameters,
  which ride charts.

## Value

A named list, empty for a penalty with no hyperparameters.

## Details

The draw is on the unconstrained scale so that whatever comes out is
admissible: a scale stays positive, a correlation stays inside its
interval, and a log-Cholesky coordinate is free on the line already. The
width is half of `sd` for the reason
[`modelterms7::term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.html)
halves it, a chart mapping a width of one onto most of its parameter's
range.

`unit` moves the center away from a one-sided bound and affects nothing
else, so a hyperparameter that is unbounded on both sides (a coordinate
of a multivariate prior's own chart) is drawn around the neutral point
of that chart whatever the term is. This is also the point that a fit
starts from; on `dr_prod` it corresponds to unit standard deviations and
no correlation.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md),
[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
