# Reject an Outer Criterion on the Expected Information Beside a Structural Term

Signals an error where `method` is an outer criterion
([`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md))
whose `information` is
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
and the model carries a structural term, and returns `NULL` invisibly
otherwise.

## Usage

``` r
structural_outer_check(method, design)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  or `NULL`.

- design:

  The design.

## Value

`NULL`, invisibly.

## Details

A criterion over a model with a structural term reads the information
over the coefficients and the term's own parameters together, and its
exact gradient and Hessian are written for the observed information of
that vector. On the expected information conditional on the past, which
[`statmod_full_information()`](https://statmodels7.github.io/statmodels7/reference/statmod_full_information.md)
computes for a filter, the determinant is a different criterion, and its
exact derivatives would need the derivative of the filter's Jacobian;
such a criterion is therefore not offered. For a term of the likelihood
shape the expected information does not exist at all (see
[`structural_expected_check()`](https://statmodels7.github.io/statmodels7/reference/structural_expected_check.md)).

## See also

[`outer_resolve()`](https://statmodels7.github.io/statmodels7/reference/outer_resolve.md),
which settles
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
as the observed information for such a model.
