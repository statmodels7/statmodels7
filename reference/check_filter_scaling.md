# Check That a Family Supports a Scaled Score

Signals an error unless the family is univariate and its second
derivative of the expected information can be evaluated at the first
observation. The derivatives of a filter driven by a scaled score read
that second derivative, and the third and fourth derivatives of the
information are finite differences of it.

## Usage

``` r
check_filter_scaling(d, y, theta)
```

## Arguments

- d:

  A distribution object.

- y:

  The response.

- theta:

  The parameters, one value per observation or one in all.

## Value

`NULL`, invisibly.
