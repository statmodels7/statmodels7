# The Coordinates a Fit Does Not Identify, Found After the Fact

The coordinates a pivoted decomposition of the penalized information at
the mode leaves out, for a fit whose own solve named none.

## Usage

``` r
deficient_coords(K)
```

## Arguments

- K:

  The penalized information at the mode, \\H + S\\.

## Value

An integer vector of coordinates of the coefficient vector, possibly
empty.

## Details

Only
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
on a pivoting decomposition reports which column it dropped, that being
a by-product of the solve that fitted the model. An `optimizers7` method
solves nothing by a pivot, and neither do the `chol`, `svd` and
`chol_crossprod` routes, so a design of less than full rank left every
coordinate unnamed and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
refused the whole matrix. The question a pivot answers is asked here
instead, on the matrix [`vcov()`](https://rdrr.io/r/stats/vcov.html)
inverts anyway, so the two cannot disagree about which model is being
reported.

The matrix is \\K = H + S\\, which is what makes a PENALIZED coordinate
safe without a clause of its own: the penalty's own curvature is in
\\S\\, so a column the design alone does not identify is identified in
\\K\\, exactly as it is identified in the augmented system the pivoted
route factorizes. Measured on two identical columns under
`ridge(~ 0 + x1 + x3)`, nothing is named.

\\K\\ is EQUILIBRATED to a unit diagonal before the pivot runs, which is
the correction
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
and the sparse rank test already carry: a smoothing parameter a
criterion sends to 1e15 and a break-point term's annealed columns both
separate the scales without flattening a direction, and per-direction
scaling forgives either. Measured, a coordinate whose curvature is 1e-14
of its neighbours' is not named.

A coordinate is a candidate only where its OWN curvature is finite and
positive, which is
[`boundary_coords()`](https://statmodels7.github.io/statmodels7/reference/boundary_coords.md)'
rule read once more: a parameter at its link's clamp makes its whole row
non-finite, and a coordinate the information carries nothing about has
an empty one. Neither is an aliasing. There the estimate stands and only
the variance does not, which is what
[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md)
handles and why the two must not be confused: aliasing reports the
coefficient itself as missing.

WHETHER THERE IS ANYTHING TO NAME IS
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)'S
VERDICT and not a tolerance of this function's, so a fit whose
information inverts is untouched and the alias can never disagree with
the variance reported beside it. That gate is load-bearing rather than
an economy. \\K\\ is \\X'X\\ up to the weights, so it SQUARES the
conditioning of the design, and a pivoted QR read at `dqrdc2`'s own
tolerance is therefore twice as strict here as it is on the augmented
system: measured on two columns made collinear to within 1e-4 – an
ordinary pair of correlated covariates – the bare pivot names one of
them while
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
inverts the matrix without difficulty, so the column would have been
thrown away for nothing. With the gate, that case is untouched and the
naming runs only where the whole matrix would otherwise have been
refused.

Which coordinate is then named is the pivot's, as it is in
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md).
Measured against that pivot where both speak, the two agree on a
duplicated column under four families and on an over-parametrized
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
term.

## See also

[`uninformative_coords()`](https://statmodels7.github.io/statmodels7/reference/uninformative_coords.md),
which holds a coordinate the information carries nothing about rather
than aliasing it, and
[`augmented_solve()`](https://statmodels7.github.io/statmodels7/reference/augmented_solve.md),
whose pivot answers the same question during the fit.
