# Draws of One Group's Effects From a Prior

Draws of One Group's Effects From a Prior

## Usage

``` r
prior_draw_groups(pen, th, D, n_draw)
```

## Arguments

- pen:

  The penalty read as a prior.

- th:

  Its hyperparameters.

- D:

  The coordinates one group carries.

- n_draw:

  The number of draws.

## Value

A matrix with `n_draw` rows and `D` columns, a group's coordinates side
by side, each row an independent draw.

## Examples

``` r
lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
pen <- penalties7::distrib_penalty(lp, n_coef = 5)
dim(statmodels7:::prior_draw_groups(pen, list(sigma = 1), 1, 12))
#> [1] 12  1
```
