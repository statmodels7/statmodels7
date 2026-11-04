# The Terms at the Vertex of Their Shape

Returns one entry per term of a fit with levels at the vertex of their
shape
([`modelterms7::term_vertex()`](https://statmodels7.github.io/modelterms7/reference/term_vertex.html)):
the equation, the term, the coefficients concerned and a sentence that
names them.

## Usage

``` r
statmod_vertex(fit)

spec_vertex(spec)
```

## Arguments

- fit:

  A
  [StatmodFit](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec:

  The specification of a fit, with its terms as the fit left them.

## Value

A list of lists with `param`, `term`, `coef` (the names as
`<param>:<coefficient>`), `part` (`"term"`, or `"interaction"` for the
interaction of
[`modelterms7::scte()`](https://statmodels7.github.io/modelterms7/reference/scte.html)
under two monotone shapes at its bound), `complete` (every level of the
term at the vertex) and `message`.

## Details

At the vertex every coordinate that the shape constrains is at its
bound, so the fitted curve is the one at the bound (constant, linear for
a plain curvature, or of slope `min_slope`) and the data go against the
shape. The columns of the block on those coordinates vanish and the
penalty is zero, so the term's smoothing parameter has no effect on the
fit.

`spec_vertex()` does the same on the terms of a fitted specification, as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
needs before the fit exists.
