# The Prior of One Random-Effect Term, for Integrating a New Group

Nodes and weights, or draws, for one group's effect under the prior the
fit estimated for it.

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

A term whose penalty is quadratic in the coefficients carries a Gaussian
prior, read off the penalty's Hessian: its first \\d\times d\\ block is
one group's precision \\\Omega\\, the block being \\I_m \otimes
\Omega\\. Any other prior, a Student t for one, is sampled with
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html),
one draw of the block giving \\m\\ independent groups.
