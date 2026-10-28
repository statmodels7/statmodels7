# Whether a Free Break-Point Has Left Its Interval

Reads the intervals that
[`pin_held_breakpoints()`](https://statmodels7.github.io/statmodels7/reference/pin_held_breakpoints.md)
recorded in the state of a design and returns `TRUE` when the position
of a free break-point lies outside its interval at the coefficients
`cf`. The objective is infinite there, so that a line search backtracks
to the inside.

## Usage

``` r
pin_infeasible(design, cf)
```

## Arguments

- design:

  A design.

- cf:

  The coefficient list, one vector for each parameter.

## Value

A single logical, `FALSE` for a design with no recorded interval.
