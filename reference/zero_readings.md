# The Optimality Check of the Coefficients a Kinked Penalty Set to Zero

The largest ratio \\\|s_j\|/\kappa_j\\ over the coefficients a penalty
with a kink set exactly to zero, where \\s_j\\ is the derivative of the
log-likelihood in that coefficient and \\\kappa_j\\ the size of the kink
there.

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

`NULL` where no coefficient is at zero, otherwise a named vector with
`ratio` and `n`, the number of such coefficients.

## Details

At a zero the penalty has no derivative, and the condition for an
optimum is not a vanishing gradient but \\\|s_j\| \le \kappa_j\\: the
likelihood's pull must not exceed what the kink holds back. A ratio at
or below one therefore says every zero is where it belongs, and a ratio
above one names a coefficient that would improve the objective by
leaving zero. No standard error enters, a coefficient at zero having
none.

A ratio of 0.89 on a lasso says the coefficient the data pull hardest is
held by the kink with 11 per cent to spare: with the other coefficients
fixed it would leave zero at a \\\lambda\\ about 11 per cent smaller,
while along a path the others move and the threshold with them. It is
not a test of significance, and values close to one are ordinary where
\\\lambda\\ was chosen on a grid, the next coefficient to enter being
close to its threshold by construction.

Measured on 200 observations of 20 columns, fitted lasso, MCP, elastic
net and standardized lasso models read 0.83 to 0.89, and a point where
one true coefficient was forced to zero and every other refitted without
it reads 20.1 while its
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
gradient is 2.7e-04.

The kink is read from the penalty's own gradient just to the right of
zero, so it is \\\lambda\\ for a lasso, \\\lambda\alpha\\ for an elastic
net and \\\lambda\\ scaled by the spread under `standardize`, with no
formula per family. The penalties are separable, so every zero is moved
off in one evaluation.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
[`certificate_readings()`](https://statmodels7.github.io/statmodels7/reference/certificate_readings.md)
