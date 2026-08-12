# The Score in a Structural Term's Own Parameters

The derivative of the weighted log-likelihood in the unconstrained
parameters of each structural term, one entry per term.

## Usage

``` r
statmod_structural_score(spec, coef, design = statmod_design(spec))
```

## Arguments

- spec:

  A
  [`StatmodSpec`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- coef:

  A named list of coefficient vectors.

- design:

  The design.

## Value

A named list of numeric vectors, one per structural term.

## Details

The filter's jacobian is the total derivative of the predictor in the
term's parameters, the recursion included, so the chain rule over the
observations is all that is needed here and no reverse pass is: what the
reverse pass answers is the other question, the derivative in the
coefficients of the equations.

## See also

[`statmod_filter_at`](https://statmodels7.github.io/statmodels7/reference/statmod_filter_at.md)
