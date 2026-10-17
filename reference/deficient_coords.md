# The Coordinates That a Fit Does Not Identify, Found After the Fact

Returns the coordinates that a pivoted decomposition of the penalized
information at the mode leaves out, for a fit whose own solve named
none.

## Usage

``` r
deficient_coords(K, held = integer(0))
```

## Arguments

- K:

  The penalized information at the mode, \\H + S\\.

- held:

  The coordinates that a kinked penalty holds at zero, from
  [`zero_kinked()`](https://statmodels7.github.io/statmodels7/reference/zero_kinked.md),
  which are not candidates.

## Value

An integer vector of coordinates of the coefficient vector, possibly
empty.

## Details

Only
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
on a pivoting decomposition reports which column it dropped, which is a
by-product of the solve that fitted the model. An `optimizers7` method
does not solve by a pivot, and the `chol` and `svd` routes do not
either, so a design of less than full rank left every coordinate unnamed
and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
signaled an error for the whole matrix. The question that a pivot
settles during the fit is settled here on the matrix that
[`vcov()`](https://rdrr.io/r/stats/vcov.html) inverts, so that the alias
and the variance refer to the same model.

The matrix is \\K = H + S\\, which makes a penalized coordinate safe
without a clause of its own: the curvature of the penalty is in \\S\\,
so a column that the design alone does not identify is identified in
\\K\\, as it is identified in the augmented system that the pivoted
route factorizes. Under `ridge(~ 0 + x1 + x3)` with two identical
columns, the function names no coordinate.

\\K\\ is equilibrated to a unit diagonal before the pivot runs, as in
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
and the sparse rank test: a smoothing parameter sent to 1e15 and the
annealed columns of a break-point term both separate the scales without
flattening a direction, and per-direction scaling accommodates both. A
coordinate whose curvature is small only because of its scale is not
named.

A coordinate is a candidate only where its own curvature is finite and
positive, which is the rule of
[`boundary_coords()`](https://statmodels7.github.io/statmodels7/reference/boundary_coords.md):
a parameter at the clamp of its link makes its whole row non-finite, and
a coordinate about which the information carries nothing has an empty
row. Both cases differ from aliasing. In them the estimate stands and
only the variance does not, a situation that
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)
handles, whereas aliasing reports the coefficient itself as missing.

[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
determines whether there is anything to name, and a tolerance of this
function does not, so a fit whose information inverts is untouched and
the alias cannot disagree with the variance reported beside it. The gate
is necessary because \\K\\ is \\X'X\\ up to the weights, so it squares
the conditioning of the design, and a pivoted QR read at the tolerance
of `dqrdc2` is twice as strict here as it is on the augmented system.
For two columns collinear to within 1e-4 (an ordinary pair of correlated
covariates), the bare pivot would name one of them while
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
inverts the matrix without difficulty. With the gate, the function
leaves such a case untouched and names a coordinate only where the whole
matrix would otherwise have been rejected.

Which coordinate is then named is determined by the pivot, as in
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md).

A coordinate that a kinked penalty holds at zero is left out of the
test, as
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
leaves it out of the matrix that it inverts. The kink identifies such a
coordinate although \\K\\ carries no curvature for it: under
`1 + lasso(~ g)` with every level of `g` coded, the direction that
raises the intercept and lowers every level is flat in \\K\\ and not in
the objective, whose minimum over it is unique. Reading the pivot on the
whole matrix would name the last level, so that coordinate is excluded.

## See also

[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md),
which holds a coordinate about which the information carries nothing
instead of aliasing it, and
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md),
whose pivot settles the same question during the fit.
