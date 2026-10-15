# The Intercept of Each Equation, on the Link Scale

The intercept-only maximum likelihood estimate, where it can be had, and
the data-based start of
[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
otherwise.

## Usage

``` r
statmod_intercepts(spec)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

## Value

A named list, one entry per distribution parameter, on the link scale;
an entry is `NULL` where neither route answered.

## Details

Two routes, tried in order.
[`distributions7::fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.html)
fits the distribution to the response with no covariates, which is the
same model with every slope set to zero and therefore a natural place
for the fit to begin; its link-scale coefficients are the intercepts.
[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
is the fallback, and its result is a list of starts, each keyed by
parameter, so a value is reached at `[[1]][[p]]`.

**The random stream is pinned and restored.** The intercept-only fit
takes several starts (`n_start = 10`) and keeps the best, and the starts
of some families are random draws. The seed is therefore fixed for the
length of this call and the caller's stream is put back afterwards, so
that the same call gives the same answer twice.

**An intercept the fit ran to the limit of its chart is replaced** by
the data-based start of
[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html).
On a log or logit chart a value past \\e^{\pm 16}\\ is a limit of the
family and not an estimate: a Student t fitted to a response no
heavier-tailed than a gaussian puts \\\nu\\ at a value far beyond any
useful estimate. A marginal criterion that estimates such a coefficient
is flat in it, and a search that starts there does not leave it, so the
data-based start is used instead. An identity chart is left alone, a
location of any size being an estimate.

## See also

[`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
