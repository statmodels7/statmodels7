# The Third Derivative of the Objective Contracted Once

\\T\[v\] = (\partial K/\partial\beta)\cdot v\\, a matrix over the
stacked coefficients.

## Usage

``` r
contract3(spec, design, d3, params, npar, offs, total, tv, key = NULL)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- d3:

  The third derivatives in the link scale, or on the expected route the
  derivative of the expected information.

- params, npar, offs, total:

  The block bookkeeping.

- tv:

  The predictors of the direction.

- key:

  A function of three parameter positions returning the name of the
  component of `d3` to read, or `NULL` for the observed route's key,
  which is symmetric in all three positions.

## Value

A square matrix.

## Details

Each block is a weighted crossproduct, the weight being \\w_i\sum_k
\ell'''\_{abk}(X_k v_k)\_i\\: the third derivative never appears as an
array.

On the expected route `d3` is the derivative of the expected information
in the predictors,
[`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html),
and the same contraction gives \\(\partial H_E/\partial\beta)\cdot v\\;
that array is symmetric in its first two positions only, so it is read
through its own `key`.
