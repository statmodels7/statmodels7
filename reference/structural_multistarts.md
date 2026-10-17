# The Starting Points of a Structural Term

Returns the starts given by
[`modelterms7::term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.html)
for the structural term of the model, when that function gives more than
one and the specification carries no values from a previous fit.

## Usage

``` r
structural_multistarts(spec, design)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

`NULL` where there is nothing to try, otherwise a list with `term`, the
unit's key, and `zeta`, the list of starts, whose first element is the
one that
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
already used.
