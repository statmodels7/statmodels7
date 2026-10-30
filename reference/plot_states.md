# Graph of the State Probabilities of a Regime Term

Selects the graph of the probability of each regime of a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term at each observation, drawn as stacked areas over the order of the
series.

## Usage

``` r
plot_states(smoothed = TRUE, legend = TRUE)
```

## Arguments

- smoothed:

  Whether to draw the smoothed probabilities (`TRUE`) or the filtered
  ones (`FALSE`).

- legend:

  Whether to draw the legend that names the elements of the graph.

## Value

A `PlotStates` object, inheriting from
[term_plot](https://statmodels7.github.io/statmodels7/reference/plot_type.md).

## Details

With `smoothed = TRUE` (the default) the probabilities are the smoothed
ones, \\P(S_t = j \mid y_1, \dots, y_n)\\, as
[`statmod_latent()`](https://statmodels7.github.io/statmodels7/reference/statmod_latent.md)
reports them. With `smoothed = FALSE` they are the filtered ones,
\\P(S_t = j \mid y_1, \dots, y_t)\\, from the forward recursion of the
term at its fitted parameters. The two coincide at the last observation
of each group. The horizontal axis is the variable named by `time` in
the term, or the row order when it has none, and a panel is drawn for
each group of a panel (at most twelve).

## See also

[`autoplot.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/autoplot.StatmodFit.md),
[`plot_regime_density()`](https://statmodels7.github.io/statmodels7/reference/plot_regime_density.md),
[`statmod_latent()`](https://statmodels7.github.io/statmodels7/reference/statmod_latent.md)

## Examples

``` r
set.seed(4)
n <- 80
s <- integer(n)
s[1] <- 1L
for (t in 2:n) s[t] <- if (runif(1) < 0.9) s[t - 1] else 3L - s[t - 1]
dd <- data.frame(t = seq_len(n), y = c(0, 3)[s] + rnorm(n, sd = 0.7))
fit <- statmod(y ~ regime(k = 2, time = t),
               distributions7::gaussian1_distrib(), dd)
ggplot2::autoplot(fit, term = "regime(k = 2, time = t)")

ggplot2::autoplot(fit, term = "regime(k = 2, time = t)",
                  plot_type = plot_states(smoothed = FALSE))
```
