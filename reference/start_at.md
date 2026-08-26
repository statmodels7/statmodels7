# Where a Fit Begins

The one generic a starting-value strategy implements. Given the model,
it returns one starting vector per distribution parameter, on the
coefficient scale, ready for the fit to begin from.

## Usage

``` r
start_at(strategy, spec, design, obj, ...)
```

## Arguments

- strategy:

  A
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md),
  which decides the method.

- spec:

  The
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md)
  being fitted.

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- obj:

  The objective, as
  [`statmod_objective()`](https://statmodels7.github.io/statmodels7/reference/statmod_objective.md)
  returns it, or `NULL`. Only
  [`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
  reads it; the other three accept `NULL`.

- ...:

  Passed to methods. No shipped method reads it.

## Value

A named list, one numeric vector per distribution parameter in the
family's order, each as long as that parameter's design is wide.

## Details

Write a method on your own subclass of
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
to add a strategy. The four shipped methods show the range:
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md)
reads only the design's widths,
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
fits a small model,
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
draws from the caller's generator, and
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
runs an optimizer over `obj`.

A method must return a full-length vector for every distribution
parameter, including the ones it has nothing to say about. Zero is the
conventional filler, and it is the value a penalized block wants anyway,
its penalty being smallest there.

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
