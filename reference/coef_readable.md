# The Quantities Reported by a Term in Place of Its Coordinates

Returns the coefficient vector of one equation, in which each term's
declared quantities replace the coordinates from which they are
computed.

## Usage

``` r
coef_readable(spec, design, fit, p, v)
```

## Arguments

- spec:

  The fitted specification.

- design:

  The design.

- fit:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- p:

  The distribution parameter naming the equation, a string.

- v:

  The named coefficient vector of that equation.

## Value

`v`, in which each term's declared quantities replace the coordinates
from which they are computed. The vector `v` is returned unchanged where
no term of the equation declares any.

## Details

A term declares its quantities through
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which returns the quantities and their Jacobian with respect to the
coefficients from which they are computed. The columns that the Jacobian
touches are exactly the coordinates to replace. A coordinate that enters
no declared quantity stays in place.

That rule keeps a developed parameter intact. A development is a vector
of coefficients over covariates with no single value to report, so the
term declares nothing for it and its coefficients stay.

The names are composed from the term's own label, as the names of its
coefficients are, so that two terms of one kind in one formula stay
apart.

## See also

[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
the caller,
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)
for what a term declares.
