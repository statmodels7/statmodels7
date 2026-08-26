# Predictors at Rows That Continue the Series

The linear predictors and the parameters at new rows for a model
carrying a structural term, with the term's own recursion carried
forward from where the fitting data left it.

## Usage

``` r
statmod_eta_continued(fit, spec, design)
```

## Arguments

- fit:

  The fitted model.

- spec:

  The specification at the new data.

- design:

  Its design.

## Value

A list shaped as
[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)'s,
without the memoized filter objects.

## Details

The static part is the ordinary assembly: each term's block reapplied at
the new rows, the offsets re-evaluated there. What cannot be reapplied
is the structural term, whose contribution at one row is the state a
recursion has reached over every row before it. The fit is therefore run
once at the observed rows to recover that state – the level and the
driving quantity at each of them, the level being the difference between
the filtered predictor and the static one – and the term is asked to
continue from there through
[`modelterms7::term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.html).

Which of the two a call asks for is decided by the response, not by the
times. New rows carrying the response are a re-reading: the filter is
run over them from the term's own seed, and that is what a caller means
by predicting a model on another series, and is why
`predict(fit, newdata = <the fitting data>)` returns the fitted values.
New rows without it are a continuation, and must come after the observed
series. A frame carrying the response on some rows only is rejected: the
two readings differ, and picking one would answer a question that was
not asked.

A term whose contribution is a likelihood mixed over latent states is
rejected: what such a term reports at an observed row is a posterior
over states, which past the data is a predictive distribution, no single
value.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
