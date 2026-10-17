# Starting Coefficients

Builds the starting coefficients of every equation: the intercepts from
the intercept-only fit of the distribution, a weighted least-squares
start for the parametric slopes of the location where it applies, and
each term's own start where the term declares one.

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

Each equation's intercept starts at the intercept-only maximum
likelihood estimate, which
[`distributions7::fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.html)
supplies: the model from which every covariate is removed is a natural
place for the model with covariates to begin, and it costs one small
fit. A block whose term declares no start of its own starts at zero,
where its penalty is smallest, and the penalized blocks are of this
kind.
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
describes the strategy in full.

[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
is the fallback for an intercept that the fit cannot supply. It returns
one list per start, and each list is keyed by parameter, so the value
for parameter `p` is `th[[1]][[p]]`.

A start that cannot be obtained is not an error: the fit still runs,
from zero on the link scale. The two routes are tried in order.

`start` is a named list of values, a
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md),
or `NULL`, which is
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).
A strategy is evaluated once, here, before the alternation between the
coefficients and the hyperparameters begins, and it must return one
numeric vector without missing values per distribution parameter, of the
length of that parameter's coefficient vector. A named list replaces the
coefficients of the parameters that it names, and may name a subset; the
other parameters start as under
[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md).
A list without names, a name that is not a parameter, and a vector of
the wrong length are errors.
