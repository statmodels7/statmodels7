# The Third Derivative of the Objective Contracted Once

\\T\[v\] = (\partial K/\partial\beta)\cdot v\\, a matrix over the
stacked coefficients.

## Usage

``` r
contract3(spec, design, d3, params, npar, offs, total, tv)
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- d3:

  The third derivatives in the link scale.

- params, npar, offs, total:

  The block bookkeeping.

- tv:

  The predictors of the direction.

## Value

A square matrix.

## Details

Each block is a weighted crossproduct, the weight being \\w_i\sum_k
\ell'''\_{abk}(X_k v_k)\_i\\: the third derivative never appears as an
array.
