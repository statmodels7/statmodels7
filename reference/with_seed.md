# Run an Expression Under a Seed

Evaluates `expr` after `set.seed(seed)` when `seed` is not `NULL`, and
restores the state of the random number generator afterwards.

## Usage

``` r
with_seed(seed, expr)
```

## Arguments

- seed:

  `NULL` or a single number.

- expr:

  The expression.

## Value

The value of `expr`.
