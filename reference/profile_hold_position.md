# Hold the Position of a Sharp Break-Point

Moves the break-point `k` of the sharp term `term` of the equation
`param` to `value` and holds the term there, the other break-points of
the term staying where they are.

## Usage

``` r
profile_hold_position(spec, param, term, k, value)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- param, term, k:

  The equation, the key of the term and the break-point.

- value:

  The position.

## Value

`spec`, with the term replaced.
