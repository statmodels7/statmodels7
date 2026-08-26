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

## Details

An entry of `par` may be a vector of the equation's own length, a single
number used for all of them, or a function of the count returning that
many values. The function form is what expresses a structured truth – a
random effect's standard deviation, a sparse vector for a lasso to find
– without a vocabulary of its own, and it is checked to have returned
the right number of finite values, since a function that answers wrongly
would otherwise be recycled into a different model.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
