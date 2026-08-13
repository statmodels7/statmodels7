# S7 Classes for the Shipped Strategies

The classes the constructors below instantiate.

## Usage

``` r
StartIntercepts(label = character(0))

StartOrigin(label = character(0))

StartRandom(
  label = character(0),
  fn = function() NULL,
  args = list(),
  center = logical(0)
)

StartSearch(
  label = character(0),
  optimizer = NULL,
  over = NULL,
  hyper = logical(0)
)
```

## Arguments

- fn:

  The generator a random start draws from.

- args:

  Further arguments to it.

- center:

  Whether the draw is added to the intercept-only start.

## Value

An S7 object inheriting from
[`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).
