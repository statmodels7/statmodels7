# Whether a Prior Is Gaussian

A penalty whose Hessian in the coefficients does not move with them and
has no kink is a Gaussian prior. The second condition matters: a Laplace
prior also has a constant Hessian, zero away from its kink.

## Usage

``` r
prior_is_gaussian(pen, th)
```

## Arguments

- pen:

  A penalty.

- th:

  Its hyperparameters.

## Value

`TRUE` or `FALSE`.

## Examples

``` r
lp <- distributions7::fixed(distributions7::laplace_distrib(), mu = 0)
statmodels7:::prior_is_gaussian(penalties7::distrib_penalty(lp, n_coef = 3),
                                list(sigma = 1))
#> [1] FALSE
```
