# The Trace of the Determinant's Movement With the Coefficients

\\u_c = \mathrm{tr}(M\\\partial K/\partial\beta_c)\\, assembled one
crossprod per distribution parameter.

## Usage

``` r
u_vector(spec, design, coef, M, params, npar, offs, total)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- coef:

  The coefficients.

- M:

  The matrix the trace is taken against.

- params:

  The distribution's parameter names.

- npar, offs, total:

  The block sizes, their offsets and the total.

## Value

A numeric vector as long as the stacked coefficients.
