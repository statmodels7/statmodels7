# Coordinates Whose Design Column Has Vanished

Returns the coordinates of the stacked coefficient vector whose column
of the design, read at the fitted coefficients, has vanished: its norm
is at most `eps` times the largest column norm of its own equation, and
the penalized information carries nothing on its diagonal either.

## Usage

``` r
vanished_coords(spec, coef, design, K)
```

## Arguments

- spec:

  The fitted specification.

- coef:

  The coefficients, one vector per distribution parameter.

- design:

  The design that
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  built.

- K:

  The penalized information at `coef`, as
  [`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
  reads it.

## Value

An integer vector of coordinates of the stacked coefficient vector,
possibly empty.

## Details

This function covers the one case that
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
leaves unnamed by design. That function reads only \\K\\, and an empty
row of \\K\\ can arise in two different ways, which the matrix does not
distinguish. At a parameter's boundary the design column is alive and
the working weight vanishes: the estimate stands and only its variance
does not, a case that
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)
handles. Where a coefficient's own column vanishes, the data do not move
with it and it is not identified at all: a break-point run out of the
data has change and position columns that are identically zero there. A
sharp step or a truncated line written into
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
can park a break-point beyond the data in this way, and without this
rule the function would name one of the two coordinates, or would name
none.

Both conditions are required. If only the design column were used, the
function would name a penalized coordinate whose column is empty (a
level of a random effect with no observations, which its prior
identifies), and if only the diagonal were used, it would name a
coordinate at a boundary. The column is compared with its own equation's
largest, since two equations' designs need not share a scale, and the
diagonal with \\\epsilon^2\\ times the largest, the square of the same
resolution.

## See also

[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
for a column that the others span,
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)
for a coordinate at a boundary.
