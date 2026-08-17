# Options for the Unpenalized Parametric Block

How the design of the linear predictor's parametric part is built: the
storage, and the contrasts for its factors.

## Usage

``` r
linpar_options(sparse = NULL, contrasts = NULL)
```

## Arguments

- sparse:

  Whether the block is a `dgCMatrix`. `NULL`, the default, leaves it to
  [`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.html),
  which settles it at build from the size of the design.

- contrasts:

  The contrasts for the block's factors, as a named list of the kind
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s
  `contrasts.arg` takes, or `NULL` for the session's
  `options("contrasts")`.

## Value

A named list, for
[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`linpar_control`. The argument and this function are named differently
on purpose: with one name for both, the argument's default would resolve
to its own promise. [`glm`](https://rdrr.io/r/stats/glm.html) and
[`glm.control`](https://rdrr.io/r/stats/glm.control.html) keep them
apart for the same reason.

## Details

It governs the IMPLICIT
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
term – the one the bare covariates of a formula collapse into, which a
caller never writes – so this is the only place its arguments can be
given. A `linpar()` written out takes them directly.

**Sparse storage.** `sparse = TRUE` builds the block through
[`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which BUILDS it sparse rather than building a dense matrix and
compressing it – the second would cost the memory the choice exists to
avoid. Measured at 20000 rows and a factor of 1000 levels, 0.002 s and
1.8 MB against
[`stats::model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s
0.100 s and 161.5 MB, the numbers identical; and a design that would be
32 GB dense builds in 0.02 s and 19 MB, which is what says there is no
dense intermediate. It pays where the formula carries a factor of many
levels and costs more than it saves on numeric covariates, whose block
is dense whatever is asked for.

**There is no rescaling here, and that is measured rather than
omitted.** Scaling the columns and carrying the coefficients back is the
remedy for a conditioning that squares, which is what forming \\X'X\\
does;
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
fits through a QR of the design and never forms it. On columns spanning
fifteen decades the raw fit and the scaled one converge in the same
number of iterations, and both agree with
[`lm`](https://rdrr.io/r/stats/lm.html) to \\10^{-14}\\. What does move
is the SCORE the fit reports, 1.5e+02 against 9.2e-05, and that is a
reading rather than an answer: the final verdict is already arbitrated
on a dimensionless scale.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.html)

## Examples

``` r
linpar_options(sparse = TRUE)
#> $sparse
#> [1] TRUE
#> 
```
