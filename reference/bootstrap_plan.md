# The Draw of One Bootstrap Replicate

Returns a function of no arguments that draws the description of one
replicate: the simulated response of a parametric replicate, or the
frequency weights of a case or a cluster replicate.

## Usage

``` r
bootstrap_plan(object, resample, spec0, design0, aside, adjust)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- resample:

  The
  [resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).

- spec0, design0:

  Its specification and design.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that a parametric draw redraws.

- adjust:

  Whether a parametric plan carries the effects drawn.

## Value

A function.

## Details

The weights multiply the weights of the fit. A cluster replicate gives
the observations of a group the number of times that the group was
drawn. The attribute `kind` of the result is `"parametric"` or
`"weights"`; `targets` lists the random-effect coordinates of a
parametric draw, and `random` says whether the effects are redrawn.
