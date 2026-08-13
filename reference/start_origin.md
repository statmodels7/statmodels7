# Start at the Origin

Every coefficient at zero on the unconstrained scale, which for each
distribution parameter is the value its link maps zero to.

## Usage

``` r
start_origin()
```

## Value

A
[`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

The name is `start_origin` and not `start_zeros` because optimizers7
already exports the second for the starting POINT of an optimizer, and
the toolkit's members share no exported name — a collision would be
reported by
[`statmodels7_conflicts`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
and, far worse, would mean that which function a user got depended on
the order the packages were attached in.

It is the plainest starting point and rarely the best: a location
parameter on the identity link then begins at zero whatever the response
is, which on a response centred at a thousand sends the run travelling.
It is here because it is the reference other strategies are judged
against, and because a model whose equations are all on a log or logit
link is not harmed by it.

## See also

[`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)

## Examples

``` r
start_origin()
#> <start> zero
```
