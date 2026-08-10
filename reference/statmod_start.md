# Starting Coefficients

An intercept-only fit per distribution parameter, with every other
coefficient at zero.

## Usage

``` r
statmod_start(spec, design, obj, start = NULL)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- obj:

  The objective.

- start:

  Optional user starting values, a named list.

## Value

A stacked numeric vector.

## Details

[`distrib_start`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
gives a starting value computed from the data rather than from the
parameter's domain, which is what makes a fit arrive in a handful of
iterations instead of spending its budget travelling. The penalized
blocks start at zero, where their penalty is smallest.
