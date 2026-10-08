# The Latent Variables of a Fitted Structural Term

The posterior summary of the latent variable a structural term of the
likelihood shape integrates over: the posterior mean and standard
deviation of each group's break-points for a marginal break-point term
([`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
with `marginal = TRUE`), and the smoothed probability of each regime at
each observation for
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html).

## Usage

``` r
statmod_latent(fit)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  whose model carries a structural term of the likelihood shape.

## Value

For a break-point term, a data frame with one row per group and
break-point: `group`, `psi`, `mean` and `sd`. For a regime term, a data
frame with one row per observation, in the order of the data, and one
column per regime, `state1`, `state2`, and so on.

## Details

The quantities come from the same decomposition the marginal likelihood
is computed on. For a break-point term it is the posterior over a
group's intervals or quadrature nodes, with the within-interval moments
those of the fitted prior truncated to it. For a regime term it is the
forward and backward recursions, which give \\P(S_t = j \mid y_1, \dots,
y_n)\\. The computation is
[`modelterms7::term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.html)'s;
this function supplies what the term cannot see, the fitted predictors
and the model's log-density.

## See also

[`modelterms7::term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.html),
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:4, each = 8), x = rep(1:8, 4))
dd$psi <- 4.5 + rep(rnorm(4, 0, 0.4), each = 8)
dd$y <- 1 + 2 * (dd$x >= dd$psi) + rnorm(32, 0, 0.3)
fit <- statmod(y ~ jump(x, psi ~ random(~1 | id), marginal = TRUE),
               distributions7::gaussian1_distrib(), dd)
statmod_latent(fit)
#>   group psi     mean         sd
#> 1     1   1 4.839775 0.10668524
#> 2     2   1 4.839775 0.10668524
#> 3     3   1 4.839775 0.10668524
#> 4     4   1 5.087046 0.07171887

# the smoothed probability of the second regime at a change of level
dr <- data.frame(t = 1:60, y = c(rnorm(30), rnorm(30, 3)))
fr <- statmod(y ~ regime(k = 2, time = t),
              distributions7::gaussian1_distrib(), dr)
round(statmod_latent(fr)$state2[c(1, 29, 30, 31, 32, 60)], 3)
#> [1] 0.00 0.00 0.00 0.09 1.00 1.00
```
