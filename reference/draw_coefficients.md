# Draw or Validate the Coefficients of a Simulation

Returns one coefficient vector per distribution parameter, drawn from a
normal where the caller gave none.

## Usage

``` r
draw_coefficients(design, params, par, sd)
```

## Arguments

- design:

  The design blocks.

- params:

  The parameter names.

- par:

  A named list, or `NULL`.

- sd:

  The standard deviation of the drawn coefficients.

## Value

A named list of numeric vectors.
