# Draw the Whole Truth of a Simulation

Every coefficient, every hyperparameter and, where the formula carries a
structural term, that term's own parameters: drawn once, then
overwritten by whatever `par` names.

## Usage

``` r
rstatmod_truth(spec, design, par, sd)
```

## Arguments

- spec:

  The specification, built against a placeholder response.

- design:

  Its design.

- par:

  A named list, or `NULL`. See
  [`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
  for what a key may be.

- sd:

  The width of the draws.

## Value

A list with `coef`, one numeric vector per distribution parameter;
`hyper`, a data frame of the hyperparameters drawn or held; and `psi`,
the structural term's own parameters or `NULL`.

## Details

The rule the whole function is written from is that **whoever knows what
a quantity means draws it**. A coefficient of a design column has no
other owner and is drawn from a normal of width `sd`. A coordinate some
penalty covers is drawn from that penalty read as a prior, so a Gaussian
random effect gives Gaussian effects and a lasso gives Laplace ones; the
prior's own scale is drawn too, on the chart its penalty carries. A
structural term's parameters are drawn by the term, which knows their
charts.

The truth is drawn once whatever `n_sim` is, so the replicates differ in
what is random and not in what is being estimated.

A hyperparameter the term holds is used as given rather than drawn:
`s(x, lambda = 2)` says what the smoothing is and the simulation says it
too.

## See also

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md),
[`penalties7::penalty_draw()`](https://statmodels7.github.io/penalties7/reference/penalty_draw.html),
[`modelterms7::term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.html)
