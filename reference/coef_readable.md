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

  The fit.

- p:

  The distribution parameter.

- v:

  The named coefficient vector of that equation.

## Value

The vector, with the quantities in place of their coordinates.

## Details

A term says what it is about through
[`term_readable`](https://statmodels7.github.io/modelterms7/reference/term_readable.html),
which gives the quantities and the Jacobian from the coefficients. The
columns that Jacobian touches are the coordinates the quantities are
read from, and they are the ones replaced; a coordinate no quantity
reads stands where it is. That is what keeps a developed parameter
intact: its development is a vector of coefficients over covariates with
no single value to report, so the term declares nothing for it and
nothing is taken away.

The names are composed as the term composes its coefficients', from its
own label, so two terms of one kind in one formula stay apart.
