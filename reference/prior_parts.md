# The Priors of the Effects a Prediction Sets Aside, by Kind

Sorts the terms `random` sets aside into Gaussian blocks, Student t
scale mixtures and other priors. See
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md).

## Usage

``` r
prior_parts(spec, design, fit, aside)
```

## Arguments

- spec, design:

  The specification and design.

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that are not `"conditional"`.

## Value

A list with `gaussian` (blocks in the shape of
[`random_blocks()`](https://statmodels7.github.io/statmodels7/reference/random_blocks.md)),
`t` (each with `members`, `Sigma` the scale matrix and `nu`) and `other`
(each with `members`, `penalty` and `theta`).

## Examples

``` r
set.seed(3)
gg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
gg$y <- 1 + gg$x + rnorm(8)[gg$g] + rnorm(40)
fit <- statmod(y ~ x + random(~ 1 | g), distributions7::gaussian1_distrib(),
               gg)
rm <- statmodels7:::random_modes(fit@spec, "zero")
lengths(statmodels7:::prior_parts(fit@spec, statmod_design(fit@spec),
                                  fit, rm))
#> gaussian        t    other 
#>        1        0        0 
```
