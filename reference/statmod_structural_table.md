# A Structural Term's Parameters, With Standard Errors

The estimated parameters of every structural term, on their own scale,
with standard errors and intervals taken from the joint observed
information.

## Usage

``` r
statmod_structural_table(fit, level = 0.95)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- level:

  The interval's level, `0.95` by default.

## Value

A data frame with one row per parameter of each structural term,
carrying the estimate, its standard error and the interval's two ends,
plus `component` and `position`: which of the term's own parameters the
quantity belongs to, and where in the parameter vector it sits, both
read off the Jacobian's support. `NULL` where the model carries no
structural term.

## Details

A structural term contributes no design block, so nothing in the
coefficient tables can report it and the values were reachable only
through `fit@structural`. The information for them exists:
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
spans the coefficients **and** the term's own parameters, and the tail
block of its inverse is their variance on the unconstrained scale.

The interval is built on that scale and mapped back through the link,
which keeps a persistence inside \\(-1, 1)\\ and a positive loading
positive.

The standard error on the parameter scale is the delta method,
\\\|h'(\zeta)\|\\\mathrm{se}(\zeta)\\. It is reported for reading, and
the interval is not built from it.

A level held because an intercept in the same equation carries it is not
estimated, and is reported with a missing standard error. Zero would
claim it had been estimated exactly.

## See also

[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
for the joint matrix inverted,
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
for the same parameters without their uncertainty
