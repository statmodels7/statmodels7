# The Prior of One Random-Effect Term, for Integrating a New Group

The description of the prior that the fit estimated for one
random-effect term: a Cholesky factor for a Gaussian prior, and the
penalty with its hyperparameters otherwise.
[`random_nodes()`](https://statmodels7.github.io/statmodels7/reference/random_nodes.md)
turns it into nodes or draws.

## Usage

``` r
random_prior(spec, design, fit, param, key)
```

## Arguments

- spec, design:

  The specification and its design.

- fit:

  The fitted model, for the hyperparameters.

- param, key:

  The term.

## Value

A list with `gaussian` (logical), `dim` and, for a Gaussian prior,
`chol` (the lower Cholesky factor of the covariance); for any other,
`pen`, `theta` and the group count `m`.

## Details

A term whose penalty is quadratic in the coefficients and has no kink
carries a Gaussian prior, read off the penalty's Hessian: its first
\\d\times d\\ block is one group's precision \\\Omega\\, the block being
\\I_m \otimes \Omega\\. For any other prior the function returns the
penalty and its hyperparameters, and
[`random_nodes()`](https://statmodels7.github.io/statmodels7/reference/random_nodes.md)
integrates over it by quadrature (one univariate prior) or by
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html),
one draw of the block giving \\m\\ independent groups. A prior that is
not Gaussian is rejected under an inverse link whose range is not
bounded. A term that shares a covariance with other terms through a
label has no prior of its own and is rejected.
