# S7 Classes for the Shipped Strategies

The five classes that
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md),
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md)
and
[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md)
instantiate. Each inherits from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
and carries its `label`; `StartRandom` adds `fn`, `args` and `center`,
`StartSearch` adds `optimizer` and `over`, and `StartFrom` adds `fit`
and `rest`. Use the constructors, which validate; these raw ones do not.

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

StartFrom(label = character(0), fit = NULL, rest = NULL)
```

## Arguments

- label:

  A single string, the name that the strategy prints.

- fn:

  The generator a random start draws from, called as `fn(n, ...)`.
  `StartRandom` only.

- args:

  A list of further arguments to `fn`. `StartRandom` only.

- center:

  A single logical: whether the draw is added to the intercept-only
  start. `StartRandom` only.

- optimizer, over:

  An optimizers7 optimizer and an optional character vector of
  coefficient labels. `StartSearch` only.

- fit, rest:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  and a
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).
  `StartFrom` only.

## Value

An S7 object inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md):
`StartIntercepts` and `StartOrigin` with the one inherited property
`label`, `StartRandom` with `fn`, `args` and `center` besides,
`StartSearch` with `optimizer` and `over`, and `StartFrom` with `fit`
and `rest`.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`start_origin()`](https://statmodels7.github.io/statmodels7/reference/start_origin.md),
[`start_random()`](https://statmodels7.github.io/statmodels7/reference/start_random.md),
[`start_search()`](https://statmodels7.github.io/statmodels7/reference/start_search.md),
[`start_from()`](https://statmodels7.github.io/statmodels7/reference/start_from.md)
for the constructors,
[`start_at()`](https://statmodels7.github.io/statmodels7/reference/start_at.md)
for the methods registered on these classes.
