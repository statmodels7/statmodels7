# The Quantities a Fit Reports, and the Map From Its Coordinates

One row per quantity the model is about, over the coefficients of every
equation and the own parameters of a structural term, with the Jacobian
from the coordinates those were estimated on.

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

This is the one place the readable view is built, so that
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
and
[`confint.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/confint.StatmodFit.md)
cannot report a quantity under one name and index it under another.

A term says what it is about through
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which gives the quantities and the Jacobian. The coordinates that
Jacobian touches are replaced by the quantities; a coordinate no
quantity reads stands where it is, with a unit row of its own. That is
what leaves a developed parameter intact: its development is a vector of
coefficients over covariates with no single value, so nothing is
declared for it and nothing is taken away.

The joint coordinate vector is the design coefficients of every equation
in order, then the **free** parameters of each structural term. Free,
because a level an intercept in the same equation carries is held and is
absent from the information the variance comes from. A quantity that
reads a held parameter is marked: its value stands, and its variance
would be that of the rest alone, so it is not reported.

## See also

[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
[`statmod_structural_table()`](https://statmodels7.github.io/statmodels7/reference/statmod_structural_table.md)
