# Is the Count of the Active Coordinates Their Rank?

`TRUE` where every equation's design has full column rank, which is
sufficient for every subset of the active coordinates to have full rank
as well, so that \\\|A\|\\ is \\\mathrm{rank}(X_A)\\ whatever the active
set turns out to be.

## Usage

``` r
design_count_exact(design)
```

## Arguments

- design:

  The design, as
  [`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
  returns it.

## Value

A single logical.

## Details

The question that
[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)
answers concerns one active set at one point, and answering it there
would cost a decomposition per evaluation. Full column rank of the whole
design answers it for every active set at once, so it is read once and
reused, which is what keeps the count free on the path it was introduced
for.

The answer is kept in the design's `eta_memo` environment, keyed on the
column counts, which is safe because that environment is built fresh
with each design. A term registering
[`modelterms7::term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.html)
recomputes its block as the coefficients move, so the same column counts
do not imply the same columns and a memo that cannot see the refresh
must not answer; there the certificate is declined outright and
[`active_rank()`](https://statmodels7.github.io/statmodels7/reference/active_rank.md)
runs instead. That route is a decomposition of \\X_A\\, where the trace
needs the cross product of the same columns.

A design carrying a structural term has no `eta_memo` either, so there
the certificate is recomputed at each call and costs one decomposition
per equation. The branch that this serves is reached only where a
penalty with a kink is present and the penalty's Hessian vanishes on the
active coordinates.

## See also

[`active_rank()`](https://statmodels7.github.io/statmodels7/reference/active_rank.md),
which answers for one active set where this declines, and
[`block_column_rank()`](https://statmodels7.github.io/statmodels7/reference/block_column_rank.md)
for the decomposition both read.
