# The Quantities That a Fit Reports, and the Map From Its Coordinates

Returns one row per quantity that the model reports, over the
coefficients of every equation and the own parameters of a structural
term, together with the Jacobian from the coordinates on which those
were estimated.

## Usage

``` r
readable_joint(spec, design, fit)
```

## Arguments

- spec:

  The fitted specification.

- design:

  The design.

- fit:

  The fit.

## Value

A list with `name`, `value`, `parameter`, `term`, `scale` (one link per
quantity), `held` (whether the quantity reads a parameter that is not
estimated), `jacobian` (quantities by joint coordinates) and `n_design`.

## Details

This function builds the readable view that
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
share.
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
builds its own readable vector from the same
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)
declarations.

A term declares its quantities through
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which returns the quantities and the Jacobian. The coordinates that the
Jacobian touches are replaced by the quantities; a coordinate that the
quantities do not read stands where it is, with a unit row of its own. A
developed parameter is therefore left intact: its development is a
vector of coefficients over covariates with no single value, so the term
declares no quantity for it and the function replaces no coordinate.

The joint coordinate vector is the design coefficients of every equation
in order, then the free parameters of each structural term. A level that
an intercept in the same equation carries is held and is absent from the
information from which the variance is computed. A quantity that reads a
held parameter is marked: its value stands, and its variance would be
only that of the remaining coordinates, so it is not reported. A
quantity that reads a coordinate past the edge of its chart is marked in
the same way
([`structural_edge_rows()`](https://statmodels7.github.io/statmodels7/reference/structural_edge_rows.md)).

## See also

[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
[`statmod_structural_table()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_table.md)
