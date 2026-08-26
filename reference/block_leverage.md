# The Per-Observation Diagonal of Each Block of a Matrix

\\G\_{ab,i} = x\_{ia}'M\_{\[a\]\[b\]}x\_{ib}\\, the diagonal of \\X_a
M\_{ab} X_b'\\ taken one observation at a time.

## Usage

``` r
block_leverage(design, M, params, npar, offs, threads = 1L)
```

## Arguments

- design:

  The design.

- M:

  The matrix the traces are taken against.

- params:

  The distribution's parameter names.

- npar, offs:

  The block sizes and their offsets.

- threads:

  How many threads the sparse route's kernel may use.

## Value

A list of lists, `G[[a]][[b]]` a vector as long as the sample.

## Details

This is the one quantity every trace against \\M\\ reduces to. For a
matrix of the design's own form, \$\$\mathrm{tr}(M\\X'WX) = \sum_i
w_i\\(XMX')\_{ii},\$\$ so a contraction of a third or fourth derivative
never has to be assembled as a \\p\times p\\ matrix in order to be
traced against \\M\\: measured at 8000 observations and 69 coefficients,
forming it and taking the trace costs 25.5 ms where the weighted sum
costs 0.031 ms, and at 20000 observations and 503 coefficients 3480 ms
against 0.066 ms.

It is computed once per evaluation point and read by the gradient's
contraction and by every pair of the Hessian.

## See also

[`u_vector()`](https://statmodels7.github.io/statmodels7/reference/u_vector.md),
[`trace_design_form()`](https://statmodels7.github.io/statmodels7/reference/trace_design_form.md)
