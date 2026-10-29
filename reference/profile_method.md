# Methods of Profiling

`profile_method` is the abstract parent of the objects that select how
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md)
obtains the profile at each value of the grid. It serves to test
membership with
[`S7::S7_inherits()`](https://rconsortium.github.io/S7/reference/S7_inherits.html).
[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md)
is the method of this version.

## Usage

``` r
profile_method(name = character(0))
```

## Arguments

- name:

  A single string, the short name of the method.

## Value

Nothing: the class is abstract and its constructor signals an error.

## See also

[`profile_grid()`](https://statmodels7.github.io/statmodels7/reference/profile_grid.md),
[`statmod_profile()`](https://statmodels7.github.io/statmodels7/reference/statmod_profile.md).

## Examples

``` r
S7::S7_inherits(profile_grid(), profile_method)
#> [1] TRUE
try(profile_method(name = "mine"))
#> Error in S7::new_object(S7::S7_object(), name = name) : 
#>   Can't construct an object from abstract class <profile_method>
```
