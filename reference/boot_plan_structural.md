# The Draw of a Replicate of a Model With a Structural Term

Returns a function of no arguments that draws a new path of the series
from the fitted model.

## Usage

``` r
boot_plan_structural(object, resample, spec0, design0)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  The
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

- spec0, design0:

  Its specification and design.

## Value

A function that returns `list(y = )`, with the attributes `kind`,
`targets`, `random` and `mode` that
[`bootstrap_engine()`](https://statmodels7.github.io/statmodels7/reference/bootstrap_engine.md)
reads.

## Details

The path comes from
[`rstatmod_eta()`](https://statmodels7.github.io/statmodels7/reference/rstatmod_eta.md),
that is from
[`modelterms7::term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.html)
run at the fitted coefficients and the fitted parameters of the term. A
score-driven term draws the response as it runs, and the level at one
time is driven by the score of the value just drawn, so the replicate is
a new realization of the process and not the fitted path with new noise.
A latent chain draws its states. The effects of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term keep their estimates.

Only the parametric scheme applies: a resample of rows or of groups
breaks the recursion of the filter at every join.
