# The Coefficients That a Marginal Criterion Estimates

The positions, in the stacked coefficient vector, of the coordinates
that
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
and
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
estimate by maximizing the criterion: in the equation of every parameter
that
[`marginal_params()`](https://statmodels7.github.io/statmodels7/reference/marginal_params.md)
names, the coefficients that are not shrunk by a penalty.

## Usage

``` r
marginal_coords(spec, design, method)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  Its design.

- method:

  An
  [`OuterMethod()`](https://statmodels7.github.io/statmodels7/reference/OuterMethod-class.md),
  or `NULL`.

## Value

A list with `where` (integer positions in the stacked vector), `param`
and `name` (the equation and the coefficient name of each position), and
`skipped`, a named character vector giving, for every parameter that has
coefficients, that the default named and that this function left out,
the reason. Every field is empty where nothing is estimated this way.

## Details

A coefficient qualifies when no penalty covers it, or when it lies in
the null space of the quadratic penalty that covers it and that null
space is spanned by coordinates. The free columns of
[`modelterms7::s()`](https://statmodels7.github.io/modelterms7/reference/s.html)
are the second case: the null space of the Demmler-Reinsch penalty is
spanned by the coordinates of the linear part.

A coefficient that a penalty with a kink covers (lasso, SCAD, MCP) does
not qualify: it stays at the joint mode, and the criterion leaves it out
of its determinant
([`laplace_pinned()`](https://statmodels7.github.io/statmodels7/reference/laplace_pinned.md)).
Beside a penalty with a kink the default names no parameter, so no
coefficient is estimated on the criterion. A parameter that the caller
names in `marginal` has its unpenalized coefficients estimated on the
criterion, which corrects a dispersion's intercept for the columns of
the mean that are fitted with it.

The following configurations are not covered yet. For each one, a
parameter that the caller names in `marginal` signals an error, and a
parameter that the default would name is left out and listed in
`skipped`, which keeps the convention of the criterion without the
argument:

- a model carrying a structural term
  ([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
  [`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)),
  whose joint fit does not hold a coefficient;

- an equation carrying a block that moves with its coefficients
  ([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
  [`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
  and the other break-point terms that are not held at fixed positions);

- a model carrying a block that is a working linearization rather than a
  Jacobian, which is the case of a sharp
  [`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
  and
  [`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html),
  in any equation: the determinant over its columns reads no curvature;

- an equation carrying a penalty whose null space is not spanned by
  coordinates, which is the case of
  [`modelterms7::te()`](https://statmodels7.github.io/modelterms7/reference/te.html).

## See also

[`marginal_params()`](https://statmodels7.github.io/statmodels7/reference/marginal_params.md),
[`statmod_hold()`](https://statmodels7.github.io/statmodels7/reference/statmod_hold.md),
which holds these coefficients in the inner fit.
