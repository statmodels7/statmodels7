# A Joint Filter's Second-Order Recursion at the Fitted Point

Runs
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
for a joint filter at the parameters at which it was fitted and returns
the Jacobian of every filtered predictor beside the second derivative
contracted against the weighted score.

## Usage

``` r
joint_filter_curvature(spec, design, f, Vs, gl, H, D3, theta)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- f:

  The filter's entry of
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)'s
  `filters`.

- Vs:

  The static rows, one matrix per distribution parameter over the
  coefficients followed by the filter's parameters.

- gl, H, D3:

  The family's first three link-scale derivatives at the fitted
  predictors.

- theta:

  The per-observation parameters at the fitted predictors.

## Value

The list that
[`modelterms7::term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.html)
returns.

## Details

The weights are the link-scale gradient of each filtered parameter times
the observation weights, so the contraction is the term of the observed
information that the curvature of the levels contributes; the driving
vector's derivatives are read at the fitted predictors.
