# The Plan of a Delete-One Refit

Returns the function of the position of a unit that gives the draw of
the refit without that unit, in the form that
[`boot_body()`](https://statmodels7.github.io/statmodels7/reference/boot_body.md)
runs.

## Usage

``` r
jackknife_plan(w0, units, refresh)
```

## Arguments

- w0:

  The prior weights of the fit.

- units:

  A list with the rows of each unit.

- refresh:

  Whether the model has a block that moves with its coefficients.

## Value

A function of the position.

## Details

The function is built here so that its environment holds only the
weights and the units, which is what a worker process receives. The
weights of the rows of the unit are set to zero. For a model with a
block that moves with its coefficients the draw carries the position as
its seed, so that the restarts of the refit do not depend on the number
of workers.
