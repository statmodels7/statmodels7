# A Structural Term's Parameters, With Standard Errors

Returns the estimated parameters of every structural term, on their own
scale, with standard errors and intervals taken from the joint observed
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

A data frame with one row per quantity that each structural term reports
and the columns `parameter`, `term`, `name`, `estimate`, `se`, `lower`,
`upper`, `held` (the quantity is a held parameter), `edge` (it depends
on a coordinate past the edge of its chart), `component` and `position`
(the one of the term's own parameters to which the quantity belongs, and
where in the parameter vector it sits; both are read off the Jacobian's
support). The value is `NULL` where the model carries no structural
term.

## Details

A structural term contributes no design block, so the coefficient tables
do not report it.
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
spans the coefficients and the term's own parameters, and the tail block
of its inverse (with the penalty Hessians added) is the variance of
those parameters on the unconstrained scale. A covariance class split
between the coefficients and the term's own parameters is placed by
[`joint_penalty_at()`](https://statmodels7.github.io/statmodels7/reference/joint_penalty_at.md),
as in
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

The quantities reported are those of
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which are not always the coordinates on which the term was estimated (a
persistence of a score-driven term is an autoregressive coefficient, a
function of the partial autocorrelations). The standard error of a
quantity is the delta method over the joint variance of the term's
parameters, \\\sqrt{J V J^\top}\\, with \\J\\ the Jacobian that
`term_readable()` returns (a persistence of a term with \\q \> 1\\ reads
several coordinates). The interval is built on the scale that the term
declares for the quantity (for a score-driven term, the log scale for a
positive loading, the scale of the partial autocorrelation's link for
the persistence when \\q = 1\\, and the identity scale for a persistence
when \\q \> 1\\): the standard error is carried to that scale by the
derivative of its link and the end points are mapped back. The interval
of a persistence is therefore confined to \\(-1, 1)\\ only when \\q =
1\\.

A level held because an intercept in the same equation carries it is not
estimated, and its standard error is reported as missing, as is that of
any quantity whose Jacobian row reads a held parameter.

A quantity that depends on a coordinate past the edge of its chart (see
[`structural_edge_rows()`](https://statmodels7.github.io/statmodels7/reference/structural_edge_rows.md))
is reported without a standard error and an interval, for the reason
that a correlation at the boundary is: the chart's derivative tends to
zero there, so the delta method reports a standard error that collapses
exactly where the quantity stops being identified.

## See also

[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
for the joint matrix inverted,
[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md)
for the same parameters without their uncertainty
