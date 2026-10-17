# The Per-Observation Diagonal of Each Block of a Matrix

Computes \\G\_{ab,i} = x\_{ia}'M\_{\[a\]\[b\]}x\_{ib}\\, the diagonal of
\\X_a M\_{ab} X_b'\\ taken one observation at a time.

## Usage

``` r
block_leverage(design, M, params, npar, offs, threads = 1L)
```

## Arguments

- design:

  The design.

- M:

  The matrix against which the traces are taken.

- params:

  The distribution's parameter names.

- npar, offs:

  The block sizes and their offsets.

- threads:

  How many threads the sparse route's kernel may use.

## Value

A list of lists, `G[[a]][[b]]` a vector as long as the sample, or `NULL`
where equation `a` or equation `b` has no coefficients.

## Details

This is the one quantity to which every trace against \\M\\ reduces. For
a matrix of the design's own form, \$\$\mathrm{tr}(M\\X'WX) = \sum_i
w_i\\(XMX')\_{ii},\$\$ so a contraction of a third or fourth derivative
never has to be assembled as a \\p\times p\\ matrix in order to be
traced against \\M\\.

It is computed once per evaluation point and read by the gradient's
contraction and by every pair of the Hessian.

Where both blocks are sparse, the product of their densities is below
\\10^{-3}\\ and the number of pairs is below \\2 \times 10^{7}\\, the
diagonal is computed by
[`leverage_pairs()`](https://statmodels7.github.io/statmodels7/reference/leverage_pairs.md)
from the nonzeros of each row. The pair count bounds the work in
absolute terms, which the density by itself does not. Otherwise the
diagonal is read off the dense product that has the fewer columns.

## See also

[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
