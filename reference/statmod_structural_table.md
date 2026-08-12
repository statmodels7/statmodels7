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
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- level:

  The interval level.

## Value

A data frame with one row per parameter of per structural term, or
`NULL` where the model carries none.

## Details

A structural term contributes no design block, so nothing in the
coefficient tables can report it and the values were reachable only
through `fit@structural`. The information for them exists:
[`statmod_full_information`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
spans the coefficients AND the term's own parameters, and the tail block
of its inverse is their variance on the unconstrained scale.

The interval is built on that scale and mapped back through the link,
which is what keeps a persistence inside \\(-1, 1)\\ and a positive
loading positive; the standard error on the parameter scale is the delta
method, \\\|h'(\zeta)\|\\\mathrm{se}(\zeta)\\, and is reported for
reading rather than for building the interval from. A level held because
an intercept in the same equation carries it is not estimated and is
reported with a missing standard error rather than a zero one.

## See also

[`statmod_full_information`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
