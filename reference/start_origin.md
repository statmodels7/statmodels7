# Start at the Origin

Every coefficient at zero on the unconstrained scale, which for each
distribution parameter is the value its link maps zero to.

## Usage

``` r
start_origin()
```

## Value

A `StartOrigin` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## When it is a poor choice

A location parameter on the identity link begins at zero whatever the
response is. On a response centered at a thousand the run then starts a
thousand units away and has to travel there, and on a badly scaled model
it may not arrive at all.
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
the default, exists for that reason.

It is harmless where every equation is on a log or a logit link, since
zero is then an ordinary interior value of the parameter, and it is the
reference the other strategies are judged against.

## The name

`start_origin`, not `start_zeros`. optimizers7 already exports
`start_zeros()` for the starting point of an optimizer, and no two
members of this toolkit export the same name. A collision would be
reported by
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md),
and, worse, which function a caller reached would depend on the order
the packages were attached in.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
the default and usually the better choice.

## Examples

``` r
start_origin()
#> <start> zero

set.seed(1)
dd <- data.frame(x = runif(40))
dd$y <- 1 + 2 * dd$x + rnorm(40, sd = 0.3)
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)

# Every coefficient of every equation is zero, intercepts included.
s <- start_at(start_origin(), spec, statmod_design(spec), NULL)
s
#> $mu
#> [1] 0 0
#> 
#> $sigma
#> [1] 0
#> 
all(unlist(s) == 0)
#> [1] TRUE
```
