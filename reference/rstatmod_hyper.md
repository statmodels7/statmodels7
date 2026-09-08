# Draw One Penalty's Hyperparameters

A value per hyperparameter, drawn on the chart the penalty carries for
it around the neutral point
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

  How far from a one-sided bound the draw is centred, passed to
  [`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md).
  Half for a prior over a structural term's own parameters, which ride
  charts.

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

`unit` moves the centre away from a ONE-SIDED bound and reaches nothing
else, so a hyperparameter unbounded on both sides – a coordinate of a
multivariate prior's own chart – is drawn around that chart's neutral
point whatever the term is. That is the point a fit starts from as well,
and on `dr_prod` it reads as unit standard deviations and no
correlation, so a simulation of a covariance block begins where a reader
would put it.

## See also

[`rstatmod_truth()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_truth.md),
[`penalty_theta_start()`](https://statmodels7.github.io/statmodels7/reference/penalty_theta_start.md)
