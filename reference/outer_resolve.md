# Settle a Marginal Criterion's Information Against the Family

Replaces `hessian = "auto"` by the information that the criterion uses
for `distrib`: `"expected"` for a Student t, a skew t or a Cauchy
response, wrapped or not, and `"observed"` for every other family. Any
other value is returned unchanged.

## Usage

``` r
outer_resolve(method, distrib)
```

## Arguments

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md).

- distrib:

  The response family.

## Value

`method`, in which `hessian` is settled.

## Details

The log-density of these families is not concave in the location, so
their observed information is not positive definite at an outlier, and
the Laplace determinant built on it can fail there. The expected
information is therefore used for these families, and the observed
information for the others.
