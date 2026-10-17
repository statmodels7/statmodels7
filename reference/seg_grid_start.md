# Choose a Break-Point Term's Starting Positions on a Grid

Runs
[`modelterms7::seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
on a
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
term whose starting positions the caller did not name, and returns the
specification unchanged for anything else.

## Usage

``` r
seg_grid_start(tm, data, response)
```

## Arguments

- tm:

  One term specification, as the formula interpreter produced it.

- data:

  The data frame against which the term is built.

- response:

  The evaluated left-hand side, or `NULL`.

## Value

`tm`, with `psi` set to the grid's choice where the rule applies, and
unchanged where it does not: a term of another kind, a term whose `psi`
the caller named, or a response against which the rule cannot score.

## Details

The objective has local optima in the break-points, and the iteration
converges from within a basin around its starting position, so the start
determines which optimum is found. The term's own default is a
conventional start, the interior quantiles of the covariate, which
depend on the covariate and not on the response. Scoring a grid of
candidate positions on the response gives a start in the basin of the
best fit of the linear part.

[`modelterms7::seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
performs one linear fit for each increasing combination of `npsi` points
of a grid of ten, `choose(10, npsi)` fits in all. The score is exact for
a gaussian response, so the rule places a starting value and does not
fit. The rule is applied in whichever equation the term sits, because
the response is the quantity against which the rule scores, even where
the term develops a scale. It is skipped where the response is not plain
numbers (a censored one, a matrix or one with missing values), because
those have no reading here.

A term whose `psi` the caller names keeps it, and naming `psi` also
turns off the grid.

## See also

[`modelterms7::seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.html),
[`statmod_terms()`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
