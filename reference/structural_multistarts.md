# The Starting Points a Structural Term Asks For

The starts
[`modelterms7::term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.html)
gives for the model's structural term, where it gives more than one and
the specification carries no values a fit arrived at.

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
unit's key, and `zeta`, the list of starts, the first being the one
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
already used.
