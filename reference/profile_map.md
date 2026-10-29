# Run Independent Units of a Profile

The one place through which a profile reaches the worker processes of
[`worker_map()`](https://statmodels7.github.io/statmodels7/reference/worker_map.md).

## Usage

``` r
profile_map(ctx, cores, n, body)
```

## Arguments

- ctx:

  The context of the call.

- cores:

  The number of worker processes.

- n:

  How many units.

- body:

  The unit's body, a function of the unit index.

## Value

The list of results, in unit order.
