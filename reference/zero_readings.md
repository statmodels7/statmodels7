# The Optimality Check of the Coefficients That a Kinked Penalty Sets to Zero

The largest ratio \\\|s_j\|/\kappa_j\\ over the coefficients that a
penalty with a kink sets exactly to zero, where \\s_j\\ is the
derivative of the log-likelihood in that coefficient and \\\kappa_j\\
the size of the kink there.

## Usage

``` r
zero_readings(spec, design, obj, coef, hyper, lab)
```

## Arguments

- spec, design, obj, coef, hyper:

  The fit's specification, design, objective, coefficients and
  hyperparameters.

- lab:

  [`coef_labels()`](https://statmodels7.github.io/statmodels7/reference/coef_labels.md)
  of the design.

## Value

`NULL` where no coefficient at zero has a finite ratio, otherwise a
named vector with `ratio` and `n`, the number of coefficients at zero.

## Details

At a zero the penalty has no derivative, and the condition for an
optimum is not a vanishing gradient but \\\|s_j\| \le \kappa_j\\: the
pull of the likelihood must not exceed what the kink holds back. A ratio
at or below one therefore means that every zero is where it belongs, and
a ratio above one identifies a coefficient that would improve the
objective by leaving zero. The standard error does not enter, because a
coefficient at zero has none.

With the other coefficients fixed and a kink proportional to
\\\lambda\\, the coefficient with the largest ratio \\r \< 1\\ leaves
zero at a \\\lambda\\ that is the fraction \\r\\ of the fitted one,
while along a path the other coefficients move and the threshold moves
with them. The ratio is not a test of significance, and values close to
one are ordinary where \\\lambda\\ was chosen on a grid, the next
coefficient to enter being close to its threshold by construction.

The kink is read from the penalty's own gradient just to the right of
zero, so it is \\\lambda\\ for a lasso, \\\lambda\alpha\\ for an elastic
net and \\\lambda\\ scaled by the spread under `standardize`, with no
formula per family. The penalties are separable, so every zero is moved
off in one evaluation.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
