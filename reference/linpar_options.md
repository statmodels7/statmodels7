# Options for the Unpenalized Parametric Block

Specifies how the parametric part of each linear predictor is built: in
what storage, and with which contrasts for its factors. Pass the result
as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`linpar_control`.

It is meant for a formula naming a factor of many levels, where
`sparse = TRUE` keeps the design in sparse storage instead of a dense
matrix with one column per level.

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

A named list with the element `sparse` (possibly `NULL`) and, when
`contrasts` is not `NULL`, the element `contrasts`, to be passed as
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)'s
`linpar_control`.

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
compressing it, so no dense intermediate is formed. Sparse storage pays
where the formula carries a factor of many levels. A numeric covariate
has no zeros to exploit, and its sparse block costs more than it saves.

## Scaling of the columns

The block is not rescaled. Scaling the columns and carrying the
coefficients back is a remedy for a conditioning that squares, as
forming \\X'X\\ does.
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
with its default `decomposition = "qr"` fits through a QR factorization
of the augmented design and does not form \\X'X\\, and its final
convergence verdict is read on a dimensionless scale, so the scale of a
column does not decide it.

## The name of the argument

The argument of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
and this function are named differently on purpose. With one name for
both, the argument's default would resolve to its own promise and R
would report *promise already under evaluation*.
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and
[`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html) are
kept apart for the same reason.

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

# Same numbers, stored in much less memory.
c(dense = class(Xd)[1], sparse = class(Xs)[1])
#>       dense      sparse 
#>    "matrix" "dgCMatrix" 
all.equal(as.matrix(Xs), Xd, check.attributes = FALSE)
#> [1] TRUE
c(dense = object.size(Xd), sparse = object.size(Xs))
#>   dense  sparse 
#> 3341296  191088 
```
