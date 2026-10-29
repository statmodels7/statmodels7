# A Specification With Named Coefficients Held

Writes the values of the labelled coefficients into `held_coef`,
replacing a hold of the same coefficient and keeping the other holds.

## Usage

``` r
profile_hold(spec, param, coefname, value)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- param, coefname:

  The equation and the name of each coefficient.

- value:

  Their values.

## Value

`spec`, with `held_coef` set.
