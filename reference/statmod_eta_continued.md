# Predictors at Rows That Continue the Series

Computes the linear predictors and the parameters at new rows for a
model that carries a structural term, whose own recursion continues from
where the fitting data left it.

## Usage

``` r
statmod_eta_continued(fit, spec, design, deriv = FALSE)
```

## Arguments

- fit:

  The fitted model.

- spec:

  The specification at the new data.

- design:

  Its design.

- deriv:

  Whether to return the forecast's derivative rows.

## Value

A list shaped as
[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)'s,
without the memoized filter objects and with `eta_static`, the static
predictor at the new rows before the continued levels are added, plus
`cont_cols`: for each equation whose filter was continued with
`deriv = TRUE`, a list with `X` (the derivative of the predictor at the
new rows in the coefficients of every equation and the term's free
parameters) and `key` (the names of its columns in the variance matrix).

## Details

The static part is the ordinary assembly: each term's block is reapplied
at the new rows, and the offsets are re-evaluated there. The structural
term cannot be reapplied, because its contribution at one row is the
state that a recursion has reached over every row before it. The fit is
therefore run once at the observed rows to recover that state, which
consists of the level and the driving quantity at each of them (the
level is the difference between the filtered predictor and the static
one), and the term continues from there through
[`modelterms7::term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.html).

Only rows without the response are continued, and they must come after
the observed series. Rows that carry it are a series of their own,
rebuilt by
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
and read as the fitting rows are, which is why
`predict(fit, newdata = <the fitting data>)` returns the fitted values.

A term whose contribution is a likelihood mixed over latent states is
rejected, because such a term reports a posterior over states at an
observed row, and past the data this posterior is a predictive
distribution and not a single value.

With `deriv = TRUE` the continuation is differentiated as well, in the
coefficients of every equation and the term's free parameters
([`continued_deriv_inputs()`](https://statmodels7.github.io/statmodels7/reference/continued_deriv_inputs.md)),
and a standard error of the forecast reads this derivative.

## See also

[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
