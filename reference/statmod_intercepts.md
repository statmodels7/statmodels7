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

## See also

[`statmod_start()`](https://statmodels7.github.io/statmodels7/reference/statmod_start.md)
