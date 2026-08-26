# The Quantities a Term Reports in Place of Its Coordinates

The coefficient vector of one equation with each term's declared
quantities put where the coordinates they are read from were.

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

`v`, with each term's declared quantities in place of the coordinates
they are read from. `v` unchanged where no term of the equation declares
any.

## Details

A term says what it is about through
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which returns the quantities and the Jacobian from the coefficients they
are read from. The columns that Jacobian touches are exactly the
coordinates to replace; a coordinate no quantity reads stays where it
is.

That rule keeps a developed parameter intact. A development is a vector
of coefficients over covariates with no single value to report, so the
term declares nothing for it and nothing is taken away.

The names are composed as the term composes its coefficients', from its
own label, so two terms of one kind in one formula stay apart.

## See also

[`coef.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/coef.StatmodFit.md),
the caller,
[`modelterms7::term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.html)
for what a term declares.
