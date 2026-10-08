# The Correction for the Estimated Hyperparameters, With Estimated Coefficients

The branch of
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
for a fit whose criterion also estimated coefficients. The correction is
the trace of the information against the part of the unconditional
variance that the hyperparameters add to the block variance of
[`marginal_vcov()`](https://statmodels7.github.io/statmodels7/reference/marginal_vcov.md).

## Usage

``` r
marginal_edf_correction(mc, coef, hyper, method, idx, expected, approx, zero)
```

## Arguments

- mc:

  The list
  [`fit_marginal_context()`](https://statmodels7.github.io/statmodels7/reference/fit_marginal_context.md)
  returns.

- coef, hyper:

  The fitted coefficients and hyperparameters.

- method:

  The
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md)
  the fit ran.

- idx:

  The hyperparameters' index.

- expected, approx:

  Passed to
  [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md).

- zero:

  The answer where the correction cannot be computed.

## Value

A list with `total`, `per` and `n_hyper`, as
[`statmod_edf_correction()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf_correction.md)
returns it.

## Details

With \\B = (B\_\eta, B\_\gamma)\\ the movement of the integrated
coefficients, \\V_v\\ the inverse of the negative Hessian of the
criterion over \\(\eta, \gamma)\\ and \\V\_\gamma\\ the one over
\\\gamma\\ alone, the correction is \$\$\mathrm{tr}\\(B V_v B^\top -
B\_\gamma V\_\gamma B\_\gamma^\top) H\\.\$\$ The part for one
hyperparameter is \\\mathrm{tr}(b_k v\_{kk} b_k^\top H)\\, with \\b_k\\
its column of \\B\_\eta\\.
