# The Degrees of Freedom That a Summary Reports

Returns the table of effective degrees of freedom that
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
reports, read on the information of the variance convention `vcov_type`,
and a note where that information is the observed one and is indefinite
at the estimate.

## Usage

``` r
summary_edf(object, vcov_type, spec, design)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- vcov_type:

  A
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object.

- spec, design:

  The specification and the design of the fit.

## Value

A list with `edf` (the table, in the shape of `fit@edf`), `info` (the
information object that the table was read on, or `NULL` where the fit
does not record it) and `note` (a string, empty where there is nothing
to report).

## Details

The table stored in the fit, `fit@edf`, was read on the information of
the outer criterion
([`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md))
or, where there is none, on the inner information; the fit records that
information in `fit@methods$edf_information`. The standard errors of a
summary read the information of `vcov_type`, so the table is recomputed
on that matrix when the two differ, and the count printed beside a
standard error is a share of the same \\H\\.

The table stored in the fit is kept in three cases: a criterion that
prices the fit with the count
([`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md))
reads the count on its own information, and a summary that recomputed it
would no longer reproduce that criterion; a model with a structural term
reads the joint vector of the coefficients and the term's parameters,
which has no expected information; and a variance convention without an
`information` property has no matrix to read.

The eigenvalues of \\(H+S)^{-1}H\\ are \\v^\top H v/(v^\top H v + v^\top
S v)\\ for the eigenvectors \\v\\ of the pencil, which lie in \\\[0,
1\]\\ for a positive semidefinite \\H\\. A direction with \\v^\top H v
\< 0\\ contributes a negative share, so the count of a term can fall
below its true value, and the total can leave \\\[0, p\]\\. The observed
information is tested for this by its Cholesky factor and, where that
fails, by the smallest eigenvalue: the matrix is indefinite when that
eigenvalue is below \\-p\\\epsilon\\\lambda\_{\max}\\, the floor with
which
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
and
[`pd_logdet_dense()`](https://statmodels7.github.io/statmodels7/reference/pd_logdet_dense.md)
separate a flat direction from rounding.

## See also

[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md),
[`vcov_information()`](https://statmodels7.github.io/statmodels7/reference/vcov_information.md)
