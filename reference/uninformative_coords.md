# Which Coordinates a Singular Information Carries Nothing About

Returns the positions whose row and column of a penalized information
are empty, so that holding them and inverting the rest reports the
variance of every other coordinate exactly.

## Usage

``` r
uninformative_coords(
  A,
  tol = 1e-10,
  share = 1 - 1e-06,
  schur = 1e-08,
  row_tol = 1e-12
)
```

## Arguments

- A:

  A penalized information, over the coefficients and any structural
  tail.

- tol:

  The relative eigenvalue below which a direction is flat.

- share:

  How much of the null space one coordinate must carry to be held, as a
  share of one.

- schur:

  The largest Schur correction, in equilibrated units, that counts as no
  disturbance to the coordinates that are kept, where the candidate has
  a positive diagonal.

- row_tol:

  The share of the matrix's largest entry below which a row is read as
  accumulated rounding, where the candidate does not.

## Value

An integer vector of positions in `A`, possibly empty.

## Details

A combination \\c'\hat\beta\\ has a finite asymptotic variance exactly
where \\c\\ is orthogonal to the null space \\\mathcal{N}\\ of \\K\\,
and then \\\mathrm{Var}(c'\hat\beta) = c'K^{+}c\\ with \\K^{+}\\ the
Moore-Penrose inverse. Reading `diag(K^+)` is therefore not enough on
its own: for a \\c\\ that is not orthogonal it returns a number where
the truth is infinite. The rule is in two parts, invert away from
\\\mathcal{N}\\ and report nothing on it, and this function is the
second part.

A candidate is a coordinate with no curvature of its own. Its diagonal
is not finite, as a parameter that has run out of its range leaves it;
or it is not positive, which is the reading of
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)
for that entry and is a statement about the coordinate regardless of its
scale; or, among the coordinates whose diagonal is positive, it carries
the whole of a null direction of the Jacobi-equilibrated matrix. The
equilibration distinguishes a flat direction from scale separation, and
requiring one coordinate to carry the whole direction excludes the
collinear case, where the null vector is \\(e_1 - e_2)/\sqrt{2}\\ and
each of the two carries half of it.

A candidate is held only where holding it disturbs nothing, and the
amount that it disturbs is the Schur correction that it removes from the
others, \\K\_{Aj}K\_{jj}^{-1}K\_{jA}\\. That test has two branches,
because it has two scales.

Where \\K\_{jj}\\ is positive the coordinate has a scale of its own and
the correction, read against \\K\_{AA}\\ in equilibrated units where
that block has a unit diagonal, is \$\$\max_k
\Bigl(\frac{\|K\_{jk}\|}{\sqrt{K\_{kk}}}\Bigr)^{2} \Big/ K\_{jj} \\\le\\
\texttt{schur}.\$\$ It is dimensionless and needs no reference: a
coordinate measured a thousand times smaller than its neighbors has a
small row and a small diagonal, and the ratio is unchanged by the
scaling. Scaling one coordinate down by \\10^{-14}\\ leaves the ratio
unchanged, and the function holds no coordinate.

Where \\K\_{jj}\\ is not positive the coordinate has no scale of its own
(the information in that direction is zero and not merely small), so
that ratio is \\0/0\\ and carries no information. The only scale left is
that of the matrix, and the question becomes whether the row sits at the
level to which its entries were accumulated: \\\max_k \|K\_{jk}\| \le
\texttt{row\\tol}\max\|K\|\\. For a count model whose dispersion has
left its range, \\K\_{jj}\\ is exactly zero and the entries of the row
of that coordinate are at the rounding level of the matrix.

Where both branches fail, dropping the coordinate would leave
\\K\_{AA}^{-1}\\, the variance conditional on that coefficient being
known, which is smaller than the marginal one; the function holds no
coordinate and the error signaled by the caller stands.

## See also

[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
which signals an error for the matrix from which these are taken, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which holds them.
