# Profile by Refitting on a Grid

Selects the exact method of
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md):
at every value of the grid the model is fitted again with the
coefficients held at that value, with the tolerances, the inner
optimizer and the outer criterion of the fit.

## Usage

``` r
ProfileGrid(name = character(0), start = character(0))

profile_grid(start = c("neighbor", "fit", "both"))
```

## Arguments

- name:

  The short name of the method, a single string. The class constructor
  `ProfileGrid()` takes it, and `profile_grid()` sets it to `"grid"`.

- start:

  One of `"neighbor"` (the default), `"fit"` and `"both"`.

## Value

A `ProfileGrid` object, inheriting from
[profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md),
with the properties `name` (`"grid"`) and `start`.

## Details

The three starts differ in the coefficients from which each refit
begins. `"neighbor"` starts the refit at a value from the refit of the
adjacent value that lies nearer to the estimate, going outward from the
estimate in both directions, so that the coefficients move little
between two refits. The first value on each side starts from the
coefficients of the fit. With several held coefficients the chains run
along the first one, one pair of chains (below and above the estimate)
for each combination of the values of the others. `"fit"` starts every
refit at the coefficients of the fit, with the held ones set to the
value of the point; the points are then independent of each other.
`"both"` runs the two and keeps, at every point, the refit with the
larger penalized log-likelihood, at twice the cost. A profile with
several local maxima, as that of a break-point has, can be reached from
one start and missed from the other.

With `cores` greater than one in
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
the chains of `"neighbor"` and the points of `"fit"` are distributed
over worker processes.

## See also

[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md),
[profile_method](https://statmodels7.github.io/statmodels7/reference/profile_method.md).

## Examples

``` r
profile_grid()
#> <statmodels7::ProfileGrid>
#>  @ name : chr "grid"
#>  @ start: chr "neighbor"
profile_grid(start = "both")
#> <statmodels7::ProfileGrid>
#>  @ name : chr "grid"
#>  @ start: chr "both"
```
