# Starting Coefficients

An intercept-only fit per distribution parameter, with every other
coefficient at zero.

## Usage

``` r
statmod_start(spec, design, obj, start = NULL)
```

## Arguments

- spec:

  The specification.

- design:

  The design.

- obj:

  The objective.

- start:

  Optional user starting values (a named list), or a
  [`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Value

A stacked numeric vector.

## Details

Each equation's intercept starts at the intercept-only MLE, which
[`distributions7::fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.html)
supplies: the model with every covariate removed is the right place for
the model with them to begin, and it costs one small fit. The penalized
blocks start at zero, where their penalty is smallest.

[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
is the fallback. It returns one list per start, each keyed by parameter,
so the value wanted is `th[[1]][[p]]`; indexing the outer list by a
parameter's name gives `NULL`, and this function did that, so every
start silently fell to zero on the link scale. On a response centered at
5.84 that put the location at 0 and sent the run traveling, which is how
a Student t fitted to iris reached a variance of \\10^7\\.

A start that cannot be obtained is not an error: the fit still runs,
from a worse place. What would be an error is not noticing, which is why
the two routes are tried in order rather than one being assumed to work.

`start` is either a named list of values, a
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md)
— which is asked once, here, before the alternation between the
coefficients and the hyperparameters begins — or `NULL` for
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
which this function did before strategies existed and still does.
