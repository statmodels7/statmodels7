# The Intercept of Each Equation, on the Link Scale

The intercept-only maximum likelihood estimate, where it can be had, and
a draw from the parameter's domain otherwise.

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
same model with every slope set to zero and therefore exactly where the
fit should begin; its link-scale coefficients are the intercepts.
[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html)
is the fallback, and its result is a list of starts, each keyed by
parameter, so a value is reached at `[[1]][[p]]`.

**The random stream is pinned and restored.** That intercept-only fit
starts from draws over the parameters' domains, so it returns a
different answer on every call where a parameter is weakly identified –
fitted to `iris`, a Student t's \\\nu\\ came back at \\e^{39}\\,
\\e^{21}\\ and \\e^{17}\\ on three consecutive runs, and
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
inherited that: the same call gave log-likelihoods of -103.49, -112.11
and -111.83. A fitting function has to give the same answer twice, so
the seed is fixed for the length of this call and the caller's stream is
put back afterwards.

Pinning makes it reproducible without making it good, one draw being one
draw; several are taken and the best kept. What would make it good is a
data-based `distrib_start` method on the univariate families, which is
the design distributions7 already documents and which only its
multivariate gaussian implements.

**An intercept the fit ran to the limit of its chart is replaced** by
the data-based start of
[`distributions7::distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.html).
On a log or logit chart a value past \\e^{\pm 16}\\ is a limit of the
family and not an estimate: a Student t fitted to a response no
heavier-tailed than a gaussian puts \\\nu\\ at \\6.1 \times 10^8\\
(`gamlss.data::film90`) or \\7.0 \times 10^{10}\\ (`abdom`). Where a
marginal criterion estimates that coefficient its search starts there,
the criterion is flat in it, and on film90 `lbfgs()` did not leave it:
100 s and `not converged`, where from the data-based \\\nu = 30\\ the
same search converges in 6.9 s to the point the observed information
reaches. An identity chart is left alone, a location of any size being
an estimate.

## See also

[`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
