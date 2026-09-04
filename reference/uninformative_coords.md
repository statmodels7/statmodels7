# Which Coordinates a Singular Information Carries Nothing About

The positions whose row and column of a penalized information are empty,
so that holding them and inverting the rest reports the variance of
every other coordinate exactly.

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
its own: for a \\c\\ that is NOT orthogonal it returns a number where
the truth is infinite. The rule is in two parts, invert away from
\\\mathcal{N}\\ and report nothing on it, and this function is the
second part.

A candidate is a coordinate with no curvature of its own. Its diagonal
is not finite, which is what a parameter run out of its range leaves
behind; or it is not positive, which is
[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md)'s
own reading of that entry and is a statement about the coordinate
whatever scale it is on; or, among the coordinates whose diagonal IS
positive, it carries the whole of a null direction of the
JACOBI-EQUILIBRATED matrix. The equilibration is what tells a flat
direction from scale separation, and requiring one coordinate to carry
the whole direction is what excludes the collinear case, where the null
vector is \\(e_1 - e_2)/\sqrt{2}\\ and each of the two carries half of
it.

**A candidate is held only where holding it disturbs nothing**, and the
amount it disturbs is the Schur correction it removes from the others,
\\K\_{Aj}K\_{jj}^{-1}K\_{jA}\\. That test has two branches, because it
has two scales.

Where \\K\_{jj}\\ is POSITIVE the coordinate has a scale of its own and
the correction, read against \\K\_{AA}\\ in equilibrated units where
that block has a unit diagonal, is \$\$\max_k
\Bigl(\frac{\|K\_{jk}\|}{\sqrt{K\_{kk}}}\Bigr)^{2} \Big/ K\_{jj} \\\le\\
\texttt{schur}.\$\$ It is dimensionless and needs no reference: a
coordinate measured a thousand times smaller than its neighbours has a
small row AND a small diagonal, and the ratio sees through both. Scaling
one down by \\10^{-14}\\ leaves the ratio at one, and nothing is held.

Where \\K\_{jj}\\ is NOT positive the coordinate has no scale of its own
– the information in that direction is zero rather than small – so that
ratio is \\0/0\\ and says nothing. The only scale left is the matrix's,
and the question becomes whether the row sits at the level its entries
were accumulated to: \\\max_k \|K\_{jk}\| \le
\texttt{row\\tol}\max\|K\|\\. Measured on a count model whose dispersion
left its range, \\K\_{jj}\\ is exactly zero, the row's largest entry is
\\1.2\times10^{-15}\\ against a matrix whose largest is 2562, and
\\\lVert K\_{j\cdot}\rVert/\lVert K\rVert\\ is \\4.4\times10^{-19}\\.

Where neither branch passes, dropping the coordinate would leave
\\K\_{AA}^{-1}\\, the variance CONDITIONAL on that coefficient being
known, which is smaller than the marginal one; nothing is held and the
caller's refusal stands.

## See also

[`solve_pd()`](https://statmodels7.github.io/statmodels7/reference/solve_pd.md),
which refuses the matrix these come from, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which holds them.
