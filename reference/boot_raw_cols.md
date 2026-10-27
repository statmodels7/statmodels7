# The Columns of the Refitted Random Effects

Returns the positions, among the stacked coefficients, of the
random-effect coordinates that a parametric replicate redraws.

## Usage

``` r
boot_raw_cols(targets, cf, params)
```

## Arguments

- targets:

  The coordinates of the draw, from `attr(draw, "targets")`.

- cf:

  The coefficients of the fit, a list by equation.

- params:

  The names of the equations, in the order of the stack.

## Value

A sorted integer vector, empty when there are no targets.
