# Settle a Fresh Start of a Structural Term Against Its Equation

The adjustments a start of a structural term's own parameters needs once
the design is known: a level held by an intercept in the same equation
starts at zero, and an unheld level that the term's own start left at
zero starts at the equation's data-based intercept. A filter also
receives a second start, kept beside the first, with each score loading
at \\0.1/\mathcal{I}\\, \\\mathcal{I}\\ the expected information of its
equation's predictor at the intercept-only fit
([`predictor_information()`](https://statmodels7.github.io/statmodels7/reference/predictor_information.md));
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fits from it where the fit from the first start does not converge.

## Usage

``` r
structural_start_fixups(spec, sst, su, fresh)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- sst:

  The design's structural state, an environment carrying `zeta` and
  `held`.

- su:

  The structural units, as `attr(design, "structural")`.

- fresh:

  A logical vector, one per unit: `TRUE` where the unit's values are a
  fresh start rather than values a fit arrived at.

## Value

`NULL`, invisibly; `sst$zeta` is modified in place.

## Details

[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
applies them to the start
[`modelterms7::term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.html)
gives, and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
applies them again to every further start
[`modelterms7::term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.html)
gives, so each start of a multistart fit is treated as the first one is.
