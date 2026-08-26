# S7 Classes for the Shipped Strategies

The three classes
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md)
and
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
instantiate. Each inherits from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
and carries its `label`; only `StartRandom` adds properties of its own.
Use the constructors, which validate; these raw ones do not.

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

StartSearch(label = character(0), optimizer = NULL, over = NULL)
```

## Arguments

- fn:

  The generator a random start draws from, called as `fn(n, ...)`.
  `StartRandom` only.

- args:

  A list of further arguments to `fn`. `StartRandom` only.

- center:

  A single logical: whether the draw is added to the intercept-only
  start. `StartRandom` only.

## Value

An S7 object inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md):
`StartIntercepts` and `StartOrigin` with the one inherited property
`label`, `StartRandom` with `fn`, `args` and `center` besides.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md)
for the constructors,
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
for the methods registered on these classes.
