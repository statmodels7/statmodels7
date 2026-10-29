# Compute a Profile

The generic through which
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
reaches a method: it receives the context built for the call and the
grid of values, and returns the list of points that
[`profile_assemble()`](https://statmodels7.github.io/statmodels7/reference/profile_assemble.md)
turns into a table.

## Usage

``` r
profile_compute(method, ctx, grid, cores)
```

## Arguments

- method:

  A
  [profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md).

- ctx:

  The context of the call, as
  [`profile_context()`](https://statmodels7.github.io/statmodels7/reference/profile_context.md)
  builds it.

- grid:

  The data frame of held values, one column per coefficient.

- cores:

  The number of worker processes.

## Value

A list with one element per row of `grid`, each a list of `par`,
`hyper`, `objective`, `loglik`, `criterion` and `converged`.
