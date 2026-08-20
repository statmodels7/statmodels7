# Search the Likelihood for a Starting Point

Runs a global search ONCE, before the fit begins, over the coordinates
where the likelihood is not convex.

## Usage

``` r
start_search(optimizer = optimizers7::sa(), over = NULL)
```

## Arguments

- optimizer:

  The optimizer to search with. Defaults to
  [`optimizers7::sa()`](https://statmodels7.github.io/optimizers7/reference/sa.html).

- over:

  Optional names of the coefficients to search over, overriding the
  choice described above.

## Value

A
[`start_strategy`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md).

## Details

**Once, and not inside the fit.** A search belongs to the starting value
and not to the scoring step. Handed to `inner_optimizer` instead — as
`chain(sa(), iwls())` — it would rerun at every hyperparameter the outer
criterion tried, which on an ordinary fit is 46 times and inside
[`cv`](https://statmodels7.github.io/statmodels7/reference/cv.md) is the
folds times the path, each time returning the same answer.

**On the likelihood alone.** The penalties are off. What a starting
value has to get right is the BASIN of the likelihood, which is the one
thing the fit will not correct by itself; the penalties enter
afterwards, when their hyperparameters are estimated, and at the probe
values they represent nobody's choice — a lasso at \\\lambda = 1\\
against an unaveraged log-likelihood empties whole blocks and would
search a different model.

**Over the coordinates where the problem is not convex**, which the
toolkit already knows how to name. They are the parameters of a
structural term
([`gas`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.html)),
the coefficients of a term that recomputes its own block
([`nl`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.html)
and the break-point terms), and each equation's intercept. Everything
else keeps the default: a smooth, a ridge or a random effect is a convex
block whose optimum the scoring step reaches from anywhere, and
searching over a thousand random-effect coefficients would spend the
whole budget on the one part of the model that does not need it. A
penalized coordinate is excluded WHEREVER IT SITS, and in particular
inside such a term: a sub-formula develops a break-point or a nonlinear
parameter over groups, and those deviations are columns of the term's
own block. They are the case the rule exists for — on the likelihood
alone nothing identifies them, the penalty being what does, so a search
over them fits each group's own points and moves AWAY from the penalized
mode. `over` overrides the choice by name.

**The hyperparameters are not searched, and cannot be from here.** The
objective is the likelihood with the penalties off, in which a
hyperparameter does not appear at all, so there is nothing for a
proposal to change. A global search over them is a search over the OUTER
criterion, where each proposal costs a full refit rather than one
likelihood evaluation, and it is already available as
`statmod(outer_optimizer = optimizers7::sa())`. Kinked penalties are
outside that too: their hyperparameter has a known upper end and is
swept by a warm-started path, which a random jump would both fail to
improve on and destroy.

## See also

[`start_intercepts`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md),
[`sa`](https://statmodels7.github.io/optimizers7/reference/sa.html),
[`chain`](https://statmodels7.github.io/optimizers7/reference/chain.html)

## Examples

``` r
start_search()
#> <start> search with simulated annealing (uniform)
start_search(optimizers7::sa(maxit = 20))
#> <start> search with simulated annealing (uniform)
```
