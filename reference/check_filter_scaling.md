# Check That a Family Supports a Scaled Score

Signals an error unless the family is univariate and carries an analytic
second derivative of its expected information, which the derivatives of
a filter driven by a scaled score reach (the fourth from a stencil on
it).

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
