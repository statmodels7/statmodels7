# The Per-Observation Pieces of the Hat Matrix of a Fit

Assembles, over the coordinates `cols` of the coefficient vector, the
stacked design of the observations, the information of each observation
and its symmetric square root.
[`leverage_blocks()`](https://statmodels7.github.io/statmodels7/reference/leverage_blocks.md)
and the corrected sandwich variances of
[`hc2_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc2_sandwich.md)
are computed from these pieces.

## Usage

``` r
leverage_parts(fit, design, cols, expected, emethod)
```

## Arguments

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- design:

  The design of the fit.

- cols:

  The positions, among the coefficients stacked by equation, of the
  coordinates that the variance covers.

- expected, emethod:

  As in
  [`statmod_information_at()`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md):
  whether the information is the expected one, and the method that
  obtains it.

## Value

A list with `params` (the names of the \\d\\ equations that have
coefficients), `n`, `d`, `Z` (a sparse matrix with \\nd\\ rows, ordered
by observation and then by equation, and one column for each element of
`cols`), `W`, `root` and `root_pinv` (arrays of dimension \\n \times d
\times d\\) and `clipped` (how many observations had a negative
eigenvalue set to zero).

## Details

Let \\d\\ be the number of distribution parameters that have at least
one coefficient. Observation \\i\\ contributes the \\d \times P\\ matrix
\\Z_i\\, whose row \\a\\ holds the row of the design of equation \\a\\
in the columns of that equation and zeros elsewhere, so that the score
of the observation in the coefficients is \\u_i = Z_i' s_i\\, with
\\s_i\\ the vector of its derivatives of the log-likelihood in the \\d\\
link-scale predictors. The matrix \\W_i\\ is the \\d \times d\\
information of the observation in those predictors (the negative second
derivative, or its expectation), times the prior weight. It follows
`expected`: the expected information (with the method of the fit's
information when that is expected, and the closed form of the family
otherwise) or the observed one. The expected information of one
observation is the covariance of its score under the model and is
positive semi-definite; the observed information of one observation can
be indefinite. Its symmetric square root \\W_i^{1/2}\\ is taken from the
eigen decomposition of the positive part (eigenvalues below \\10^{-12}\\
times the largest in absolute value are set to zero), and the number of
observations with a negative eigenvalue is returned as `clipped`. The
pseudo-inverse of the root is kept beside it, since an observation can
have a singular \\W_i\\.
