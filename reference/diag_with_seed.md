# Evaluate an Expression under a Seed

Evaluates `expr` after `set.seed(seed)` when `seed` is not `NULL`, and
restores the random stream of the caller afterwards.

## Usage

``` r
diag_with_seed(seed, expr)
```

## Arguments

- seed:

  A single number or `NULL`.

- expr:

  The expression.

## Value

The value of `expr`.
