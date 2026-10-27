# A Newton Step Without the Penalty of a Term

A Newton Step Without the Penalty of a Term

## Usage

``` r
unpenalized_estimates(sys)
```

## Arguments

- sys:

  The list from
  [`penalty_system()`](https://statmodels7.github.io/statmodels7/reference/penalty_system.md).

## Value

The unpenalized estimates of the drawn coefficients, `NA` where the
coefficient is out of the step.
