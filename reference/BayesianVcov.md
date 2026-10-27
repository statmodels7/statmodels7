# S7 Classes for the Shipped Variance Conventions

These are the classes that
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md),
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md)
and
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
instantiate. Each inherits from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
and carries its `name`. Use the constructors, which set the property.

## Usage

``` r
BayesianVcov(name = character(0), information = NULL)

FrequentistVcov(name = character(0), information = NULL, sandwich = NULL)

UnconditionalVcov(name = character(0), information = NULL)

BootstrapVcov(name = character(0), n = NULL, resample = NULL, refit = NULL)
```

## Arguments

- name:

  A single string, the short name of the convention.

- information:

  `NULL` or an object that inherits from
  [distributions7::Information](https://statmodels7.github.io/distributions7/reference/Information.html):
  [`distributions7::observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.html)
  or
  [`distributions7::expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.html).
  `NULL` is resolved for each fit by the family (see
  [`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)).

## Value

An S7 object inheriting from
[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md),
with the properties `name` and `information`.

## See also

[vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
for the comparison of the three conventions.

## Examples

``` r
S7::S7_inherits(bayesian_vcov(), vcov_method)
#> [1] TRUE
class(bayesian_vcov())
#> [1] "statmodels7::BayesianVcov" "statmodels7::vcov_method" 
#> [3] "S7_object"                
bayesian_vcov()@information
#> NULL
```
