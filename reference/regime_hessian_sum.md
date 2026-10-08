# The Second-Derivative Components of a Mixture, Averaged Over Its States

For a likelihood mixed over states, each second-derivative component of
the family averaged over the states with their smoothed probabilities,
\\\bar h\_{ab,i} = \sum_k \gamma\_{ik} h\_{ab}(\theta\_{ik})\\.

## Usage

``` r
regime_hessian_sum(spec, ev, expected, approx)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- ev:

  What
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns, carrying `regimes`.

- expected, approx:

  As in
  [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md).

## Value

A named list of numeric vectors of length `n_obs`, keyed as
[`distributions7::distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.html)
keys its result.

## Details

The complete-data information is \\\sum_k X_a' \mathrm{diag}(w \gamma_k
h\_{ab,k}) X_b = X_a' \mathrm{diag}(w \bar h\_{ab}) X_b\\, so one cross
product per block replaces one per state and block. The states are the
latent components of a structural term: the regimes of
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
and the quadrature nodes or side patterns of a marginal break-point
term.
