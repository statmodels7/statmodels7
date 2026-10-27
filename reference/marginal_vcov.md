# The Variance Where the Criterion Estimates Coefficients

Returns the variance of the coefficients of a fit whose marginal
criterion estimated some of them
([`marginal_coords()`](https://statmodels7.github.io/statmodels7/reference/marginal_coords.md)),
in the block form of the penalized information with the curvature of the
criterion in place of the Schur complement of those coefficients.

## Usage

``` r
marginal_vcov(object, design, A, keep, convention)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  Its design.

- A:

  The penalized information over the kept coordinates.

- keep:

  The kept coordinates, a logical vector over the stacked coefficients.

- convention:

  The `name` of a
  [vcov_method](https://statmodels7.github.io/statmodels7/reference/vcov_method.md)
  object, `"bayesian"` or `"unconditional"`.

## Value

The variance over the kept coordinates, or `NULL` where the fit
estimated no coefficient on its criterion, where the model carries a
structural term, or where a piece cannot be read.

## Details

Write \\\gamma\\ for the estimated coefficients, \\u\\ for the rest and
\\K\\ for the penalized information. The coefficients \\u\\ are the mode
of the penalized likelihood at \\\gamma\\, and \\b_j = e_j -
K\_{uu}^{-1}K\_{u\gamma}e_j\\ is how that mode moves with \\\gamma_j\\.
With \\v\\ the outer vector, \$\$\mathrm{Var}(\hat\beta) =
K\_{uu}^{-1} + B\\V_v B^\top,\qquad V_v =
\big\[-\nabla^2\_{vv}\ell_M\big\]^{-1},\$\$ where \\K\_{uu}^{-1}\\ is
padded with zeros on \\\gamma\\ and \\B\\ holds one direction per
coordinate of \\v\\. For
[`bayesian_vcov()`](https://statmodels7.github.io/statmodels7/reference/bayesian_vcov.md)
\\v = \gamma\\ and the hyperparameters are held; for
[`unconditional_vcov()`](https://statmodels7.github.io/statmodels7/reference/unconditional_vcov.md)
\\v = (\eta, \gamma)\\ and the hyperparameters' directions are
\\-K\_{uu}^{-1}\partial^2\rho/\partial u\\\partial\eta\\.

Replacing \\-\nabla^2\_{\gamma\gamma}\ell_M\\ by the Schur complement
\\K\_{\gamma\gamma} - K\_{\gamma u}K\_{uu}^{-1}K\_{u\gamma}\\ gives back
\\K^{-1}\\, which is the convention of a fit that reads every
coefficient at the joint mode.
