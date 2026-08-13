# Where a Fit Begins

The one generic a starting-value strategy implements: given the model,
it returns one starting vector per distribution parameter.

## Usage

``` r
start_at(strategy, spec, design, obj, ...)
```

## Arguments

- strategy:

  A
  [`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

- spec:

  The specification.

- design:

  The design.

- obj:

  The objective.

- ...:

  Passed to methods.

## Value

A named list, one numeric vector per distribution parameter, each as
long as that parameter's design.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + 2 * dd$x + rnorm(40, sd = 0.3)
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
design <- statmod_design(spec)
# a strategy that needs no objective is asked with NULL for it; the ones
# that search the likelihood are handed the objective by statmod()
start_at(start_origin(), spec, design, NULL)
#> $mu
#> [1] 0 0
#> 
#> $sigma
#> [1] 0
#> 
start_at(start_intercepts(), spec, design, NULL)
#> $mu
#> [1] 2.048898 0.000000
#> 
#> $sigma
#> [1] -0.560959
#> 
```
