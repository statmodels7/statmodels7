# Options for the Unpenalized Parametric Block

Says how the parametric part of each linear predictor is built: in what
storage, and with which contrasts for its factors. Pass the result as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`linpar_control`.

The one case it exists for is a formula naming a factor of many levels,
where `sparse = TRUE` turns a design that would be gigabytes into one
that is megabytes.

## Usage

``` r
linpar_options(sparse = NULL, contrasts = NULL)
```

## Arguments

- sparse:

  A single logical, or `NULL`. `TRUE` makes the block a `dgCMatrix`,
  `FALSE` a base matrix. `NULL`, the default, leaves it to
  [`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html),
  which settles it at build time from the size of the design.

- contrasts:

  The contrasts for the block's factors, a named list of the kind
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)'s
  `contrasts.arg` takes, or `NULL` for the session's
  `options("contrasts")`. Carried on the specification, so a fold of
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  reproduces them instead of re-reading the option.

## Value

A named list with elements `sparse` and `contrasts`, to be passed as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`linpar_control`.

The argument and this function are named differently on purpose. With
one name for both, the argument's default would resolve to its own
promise and R would report *promise already under evaluation*.
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and
[`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html) are
kept apart for the same reason.

## Which term it reaches

The **implicit**
[`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
term: the one the bare covariates of a formula collapse into, which a
caller never writes and so has no other way to configure. A `linpar()`
written out in the formula takes these arguments directly and ignores
this.

## Sparse storage

`sparse = TRUE` builds the block through
[`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html),
which **builds** it sparse instead of building a dense matrix and
compressing it. The second would cost the memory the choice exists to
avoid.

Measured at 20000 rows and a factor of 1000 levels, 0.002 s and 1.8 MB
against
[`stats::model.matrix`](https://rdrr.io/r/stats/model.matrix.html)'s
0.100 s and 161.5 MB, the numbers identical; and a design that would be
32 GB dense builds in 0.02 s and 19 MB, which settles that there is no
dense intermediate. It pays where the formula carries a factor of many
levels and costs more than it saves on numeric covariates, whose block
is dense whatever is asked for.

**There is no rescaling here, and that is measured rather than
omitted.** Scaling the columns and carrying the coefficients back is the
remedy for a conditioning that squares, as forming \\X'X\\ does;
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
fits through a QR of the design and never forms it. On columns spanning
fifteen decades the raw fit and the scaled one converge in the same
number of iterations, and both agree with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) to \\10^{-14}\\. What
does move is the score the fit reports, 1.5e+02 against 9.2e-05, and
that is a reading rather than an answer: the final verdict is already
arbitrated on a dimensionless scale.

## See also

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for where it is passed,
[`modelterms7::linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.html)
for the term it configures.

## Examples

``` r
linpar_options(sparse = TRUE)
#> $sparse
#> [1] TRUE
#> 

# The whole point: a factor of many levels.
set.seed(1)
n <- 2000
dd <- data.frame(y = rnorm(n), g = factor(sample(200, n, replace = TRUE)))
dense  <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), dd,
                       linpar = linpar_options(sparse = FALSE))
sparse <- statmod_spec(y ~ g, distributions7::gaussian1_distrib(), dd,
                       linpar = linpar_options(sparse = TRUE))
Xd <- statmod_design(dense)$mu$X
Xs <- statmod_design(sparse)$mu$X

# Same numbers, one two orders of magnitude smaller.
c(dense = class(Xd)[1], sparse = class(Xs)[1])
#>       dense      sparse 
#>    "matrix" "dgCMatrix" 
all.equal(as.matrix(Xs), Xd, check.attributes = FALSE)
#> [1] TRUE
c(dense = object.size(Xd), sparse = object.size(Xs))
#>   dense  sparse 
#> 3341296  191088 
```
