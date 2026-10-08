# Coordinates at the Edge of Their Chart That Point Inward

Checks the first-order condition for a maximum at every free coordinate
that has run past the edge of a chart mapping onto a bounded set,
reading it on the bounded scale rather than on the free one. Returns the
coordinates where moving back inside raises the log-likelihood, so the
point the fit stopped at is not a maximum.

## Usage

``` r
edge_violations(spec, design, coef, edge = 8, tol = mode_error_limit())
```

## Arguments

- spec, design, coef:

  The specification, its design and the coefficients by parameter, as
  [`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  holds them at the end of a fit.

- edge:

  The free value past which a coordinate is at the edge of its chart,
  the same default as
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)'s.

- tol:

  The gain, in log-likelihood units, above which the point is not a
  maximum:
  [`mode_error_limit()`](https://statmodels7.github.io/statmodels7/reference/mode_error_limit.md).

## Value

A data frame with one row per violating coordinate and columns `kind`
(`"structural"` or `"coefficient"`), `param`, `term`, `name`, `eta`,
`target` and `gain`; zero rows where every coordinate at an edge points
outward. The design's structural state is left as it was found.

## Details

At a coordinate \\\eta\\ whose chart \\\theta = h(\eta)\\ saturates,
\\\partial\ell/\partial\eta = (\partial\ell/\partial\theta)\\h'(\eta)\\,
and \\h'\\ tends to zero at the edge whatever
\\\partial\ell/\partial\theta\\ is. A fit can therefore report a
vanishing score, and a certificate a vanishing mode error, at a point
where the log-likelihood still rises towards the interior. The condition
for a maximum on the boundary of a bounded parameter (Karush, Kuhn and
Tucker) is that the derivative on the bounded scale points outward, and
that is what is checked.

The derivative is read as a one-sided difference: each coordinate beyond
`edge` in absolute value is moved to \\\pm\\`edge`, the rest held, and
the log-likelihood is evaluated there. A coordinate is reported where
the gain exceeds `tol`, in log-likelihood units. The coordinates checked
are those
[`modelterms7::term_charted()`](https://statmodels7.github.io/modelterms7/reference/term_charted.html)
names for a structural term and for
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
and the intercept of an equation that has no other column and carries a
link that is not the identity. None of them carries a penalty, so the
log-likelihood difference is the objective's.

## See also

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
which reports a violation.
