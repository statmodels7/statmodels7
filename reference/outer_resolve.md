# Settle a Marginal Criterion's Information Against the Family

Replaces
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
by the information that the criterion uses for `distrib`:
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
for a Student t, a skew t or a Cauchy response, wrapped or not, and
[`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
for every other family. The approximation that the
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
object carries passes to the
[`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html)
object. Any other information is returned unchanged.

## Usage

``` r
outer_resolve(method, distrib, design = NULL)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- distrib:

  The response family.

- design:

  The design, or `NULL`.

## Value

`method`, in which `hessian` is settled.

## Details

The log-density of these families is not concave in the location, so
their observed information is not positive definite at an outlier, and
the Laplace determinant built on it can fail there. The expected
information is therefore used for these families, and the observed
information for the others.

For a model with a structural term the observed information is used for
every family, since the criterion over the joint vector is written on it
(see
[`structural_outer_check()`](https://statmodels7.github.io/statmodels7/reference/structural_outer_check.md)).
