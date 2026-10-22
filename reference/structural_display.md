# The Equation and the Reported Name of a Structural Term's Parameters

Returns, for names of a structural term's parameters or of its reported
quantities, the equation to which each belongs and the name under which
it is reported: the term's label, a dot, and the name.

## Usage

``` r
structural_display(tm, param, nm)
```

## Arguments

- tm:

  The structural term.

- param:

  The equation in which the term is kept.

- nm:

  The names.

## Value

A list with `param` and `name`, character vectors like `nm`.

## Details

A term of one equation belongs to that equation. A joint filter spans
several, and its names carry the parameter of their term as a prefix
(`sigma.kappa1`): the prefix gives the equation and is replaced by the
label of that term, so the loading of the scale's filter is reported as
`gas.kappa1` among the coefficients of `sigma`, as it would be without
the joint filter.
