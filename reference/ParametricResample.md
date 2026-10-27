# S7 Classes for the Shipped Resampling Schemes

These are the classes that
[`parametric_resample()`](https://statmodels7.github.io/statmodels7/reference/parametric_resample.md),
[`case_resample()`](https://statmodels7.github.io/statmodels7/reference/case_resample.md)
and
[`cluster_resample()`](https://statmodels7.github.io/statmodels7/reference/cluster_resample.md)
instantiate. Each inherits from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md).
Use the constructors, which set the properties.

## Usage

``` r
ParametricResample(name = character(0))

CaseResample(name = character(0))

ClusterResample(name = character(0), cluster = NULL, data = NULL)
```

## Arguments

- name:

  A single string, the short name of the scheme.

- cluster:

  For a cluster scheme, the one-sided formula that names the variable
  that identifies the groups.

- data:

  For a cluster scheme, the data frame that holds that variable, or
  `NULL`.

## Value

An S7 object inheriting from
[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md),
with the property `name`, and `cluster` and `data` for a cluster scheme.

## See also

[resample_method](https://statmodels7.github.io/statmodels7/reference/resample_method.md)
for the comparison of the schemes.

## Examples

``` r
class(case_resample())
#> [1] "statmodels7::CaseResample"    "statmodels7::resample_method"
#> [3] "S7_object"                   
cluster_resample(~ g)@cluster
#> ~g
#> <environment: base>
```
