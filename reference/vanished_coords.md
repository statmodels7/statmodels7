# Coordinates Whose Design Column Has Vanished

The coordinates of the stacked coefficient vector whose column of the
design, read at the fitted coefficients, has vanished: its norm is at
most `eps` times the largest column norm of its own equation, and the
penalized information carries nothing on its diagonal either.

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

  The design
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

It answers the one case
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
leaves unnamed by design. That function reads \\K\\ alone, and an empty
row of \\K\\ is two different things a matrix cannot tell apart. At a
parameter's boundary the design column is alive and the working weight
vanishes: the estimate stands and only its variance does not, which is
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)'
business. Where a coefficient's own column vanishes, nothing in the data
moves with it and it is not identified at all – a break-point run out of
the data, whose change and position columns are identically zero there.
Measured on a sharp step and a truncated line written into
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
the fit parked the break-point at 1.9e+12 with the change's column at
2e-67 and the position's at 2e-23, and without this rule one of the two,
or neither, was named.

Both conditions are required. The design column alone would name a
penalized coordinate whose column is empty – a level of a random effect
with no observations, which its prior identifies – and the diagonal
alone would name a coordinate at a boundary. The column is compared with
its own equation's largest, since two equations' designs need not share
a scale, and the diagonal with \\\epsilon^2\\ times the largest, the
square of the same resolution.

## See also

[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
for a column the others span,
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)
for a coordinate at a boundary.
