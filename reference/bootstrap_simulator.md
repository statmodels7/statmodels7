# A Simulator of the Response From a Fitted Model

Returns a function of no arguments that draws one response at the
fitting rows: the effects of the
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
terms in `aside` from the prior the fit estimated, every other
coefficient at its estimate, and the response from the family at the
parameters that gives.

## Usage

``` r
bootstrap_simulator(object, spec0, design0, aside)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec0, design0:

  Its specification and design at the fitting rows.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  whose effects are redrawn, the others keeping their estimates.

## Value

A function returning a numeric vector of `spec0@n_obs` values.

## Details

The effects are drawn with
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html),
the penalty read as a prior at the fit's hyperparameters, and placed by
the map
[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
uses, so a structured covariance (an AR(1), a compound symmetry), a
covariance a label shares between terms and a prior that is not Gaussian
(a multivariate Student t, a Laplace) are drawn alike. A smooth or a
ridge is not a random effect and keeps its estimate.

## Examples

``` r
set.seed(3)
gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
               gg)
rm <- statmodels7:::random_modes(fit@spec, "zero")
draw <- statmodels7:::bootstrap_simulator(fit, fit@spec,
                                          statmod_design(fit@spec), rm)
length(draw())
#> [1] 40
```
