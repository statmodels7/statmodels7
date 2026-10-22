# The Derivatives of the Scaling Matrix in the Predictors

Returns a function of a multiset of predictor indices that gives the
derivative of \\S = \mathcal{I}\_{GG}^{-d}\\ in those predictors, one
matrix per observation, for \\d\\ equal to 1/2 or 1.

## Usage

``` r
joint_scale_derivs(spec, theta, G, d, order)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- theta:

  The per-observation parameters at the fitted predictors.

- G:

  The positions of the filtered parameters.

- d:

  The power, 1/2 or 1.

- order:

  The highest order of derivative of \\\mathcal{I}\\ read.

## Value

A function of an integer vector of predictor positions.

## Details

With \\M = \mathcal{I}\\ for \\d = 1\\ and \\M = \mathcal{I}^{1/2}\\ for
\\d = 1/2\\, \\S = M^{-1}\\, and Leibniz's rule over the positions of
the multiset \\B\\ on \\MS = 1\\ gives \$\$S_B = -S \sum\_{\emptyset \ne
A \subseteq B} M_A S\_{B \setminus A}.\$\$ For the symmetric root the
same rule on \\XX = \mathcal{I}\\ gives \$\$X X_B + X_B X =
\mathcal{I}\_B - \sum\_{\emptyset \ne A \subsetneq B} X_A X\_{B
\setminus A},\$\$ a Sylvester equation, which is solved in the
eigenbasis of \\\mathcal{I}\\ by dividing by \\\sqrt{\lambda_i} +
\sqrt{\lambda_j}\\. The derivatives of \\\mathcal{I}\\ are the family's
([`distributions7::distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.html)
and its siblings, on the link scale). Every derivative is memoized by
its multiset.
