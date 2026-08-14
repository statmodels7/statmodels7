# Choose a Break-Point Term's Starting Positions on a Grid

Runs
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.html)
on a
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
term whose starting positions the caller did not name, and returns the
specification unchanged for anything else.

## Usage

``` r
seg_grid_start(tm, data, response)
```

## Arguments

- tm:

  One term specification.

- data:

  The data frame the term is built against.

- response:

  The evaluated left-hand side, or `NULL`.

## Value

The specification, with `psi` set where the rule applies.

## Details

The objective has local optima in the break-points and the iteration
converges from within a basin around where it starts, so where a run
begins decides what it finds. Measured over eight samples and four
starting positions on a joint jump and change of slope, the fraction of
runs recovering the break-point is 0 to 0.5 from a single conventional
start and 1 from the grid. The term's own default is a conventional
start: the interior quantiles of the covariate, which look at the
covariate and not at the response.

The rule costs `k` linear fits and is exact for a gaussian response, so
it places a starting value and does not fit. Two things are therefore
not asked of it. It is applied whatever equation the term sits in, the
response being what there is to score against even where the term
develops a scale; and it is skipped where the response is not plain
numbers – a censored one, or a matrix – rather than being given a
reading of its own.

A caller who names `psi` has said where to begin and is left alone,
which is also how the grid is turned off.

## See also

[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.html),
[`statmod_terms`](https://statmodels7.github.io/statmodels7/reference/statmod_terms.md)
