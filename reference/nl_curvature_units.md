# The nl() Blocks Whose Full Curvature an Inner Fit Uses

Returns the refreshable units of the design that belong to an
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
term whose `curvature` setting is `"full"` or `"auto"`, each with that
setting as the element `curvature`.

## Usage

``` r
nl_curvature_units(design, all = FALSE)
```

## Arguments

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

- all:

  `TRUE` to include the units whose setting is `"gauss_newton"` as well.

## Value

A list of units, possibly empty.

## Details

A frozen unit (a sharp break-point's working linearization) has no
second derivative and is not among them, and neither is a unit of
another term. A term built before the setting existed is read as
`"auto"`.
