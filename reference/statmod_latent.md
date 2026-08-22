# The Posterior Break-Points of a Marginal Term

The posterior mean and standard deviation of each group's latent
break-points in a fitted model carrying a marginal break-point term
([`jump`](https://statmodels7.github.io/modelterms7/reference/jump.html),
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.html) or
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
with `marginal = TRUE`).

## Usage

``` r
statmod_latent(fit)
```

## Arguments

- fit:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  whose model carries a structural term of the likelihood shape.

## Value

A data frame with one row per group and break-point: `group`, `psi`,
`mean` and `sd`.

## Details

The quantities come from the same decomposition the marginal likelihood
is computed on: the posterior over a group's intervals or quadrature
nodes, with the within-interval moments those of the fitted prior
truncated to it. The computation is
[`term_latent`](https://statmodels7.github.io/modelterms7/reference/term_latent.html)'s;
this function supplies what the term cannot see, the fitted predictors
and the model's log-density.

## See also

[`term_latent`](https://statmodels7.github.io/modelterms7/reference/term_latent.html),
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)

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
#> 1     1   1 4.839823 0.10665307
#> 2     2   1 4.839823 0.10665307
#> 3     3   1 4.839823 0.10665307
#> 4     4   1 5.087019 0.07169724
```
