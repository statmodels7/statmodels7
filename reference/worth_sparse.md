# Is a Matrix Worth Factorizing Sparsely?

Whether a symmetric matrix is large enough and sparse enough that a
sparse Cholesky beats a dense one, asked of the MATRIX and of nothing
else.

## Usage

``` r
worth_sparse(M, min_dim = 100L, max_density = 0.1)
```

## Arguments

- M:

  A symmetric matrix.

- min_dim:

  The smallest order worth the fixed cost.

- max_density:

  The largest fraction of nonzeros worth it.

## Value

A single logical.

## Details

The two conditions are the measured crossover and not a preference. On
the penalized information of a random intercept over \\m\\ levels at
20000 observations, the sparse route against the dense one – coercion,
factorization, log-determinant and full inverse, each timed with the
repetition loop sized by elapsed time:

|       |       |             |                 |             |
|-------|-------|-------------|-----------------|-------------|
| **m** | **p** | **density** | **whole route** | **inverse** |
| 20    | 23    | 0.282       | 1.08x           | **0.13x**   |
| 50    | 53    | 0.128       | 1.06x           | **0.33x**   |
| 100   | 103   | 0.067       | 1.18x           | 1.14x       |
| 200   | 203   | 0.034       | 1.74x           | 2.9x        |
| 500   | 503   | 0.014       | 5.28x           | 7.6x        |
| 1000  | 1003  | 0.007       | 2.50x           | 11.2x       |
| 2000  | 2003  | 0.003       | 4.20x           | 7.0x        |

The whole-route column builds the matrix afresh on each repetition and
is the one the thresholds are read from. The inverse column is the
like-for- like comparison, each route carrying its own factorization; an
earlier version of it timed the sparse solves against a factor built
once outside the loop, which flattered the sparse side without changing
where it loses.

Below about a hundred coefficients the sparse route LOSES, and loses
badly: its fixed cost is the coercion and the S4 dispatch around it,
which does not shrink with the matrix. On the fully dense penalized
information of a single smooth (p = 16, density 1) it measures 0.01x, a
hundred times slower, which is what the size condition is there to
prevent.

**Both quantities are read off the matrix, and the first one is its
STORAGE.** A matrix held as a base matrix is refused whatever its zeros,
which reads like a test of the container rather than of the mathematics,
so it is worth saying why it is neither an oversight nor a term test.
[`statmod_information_at`](https://statmodels7.github.io/statmodels7/reference/statmod_information_at.md)
accumulates into the design's own kind, so the penalized matrix is
stored sparsely exactly when the design is, and modelterms7 builds a
block sparse only when asked (`sparse = TRUE`, whose default is
`FALSE`). Measured on `y ~ 0 + g + s(x)` over 400 levels at 20000
observations, whose penalized matrix is 5 per cent nonzero either way:
built dense the fit takes 104.24 s and this factorization is **0.16 per
cent** of it, the time being in the \\O(np^2)\\ products against a dense
design (`statmod_information_at` 48.8 per cent, `crossprod` 57.2 per
cent of self time); built sparse the same fit takes 2.19 s. So where the
storage is dense the factorization is not what a fit is spending its
time on, and coercing a dense \\p \times p\\ matrix here to save a share
of that size would cost more than it returns.

**The like-for-like comparison is the one that says this is not a term
test**, and it is the check `piano_lme4.txt` section 5 asks for. With
every design built the same way, this route is worth 1.38x on
`0 + g + s(x)` over 400 levels, 1.33x on `random(~1|g)` over 500 and
1.07x on `s(x, by = g)` over 60 – an unpenalized indicator block, a
random effect and a factor-`by` smooth, gaining together and in the
order their sizes predict. Nothing here asks which term or which family
produced the matrix.

## See also

[`pd_factor`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
