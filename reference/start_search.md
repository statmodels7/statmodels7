# Search the Likelihood for a Starting Point

Runs a global search once, before the fit begins, over the coordinates
where the likelihood is not convex. Built for a model whose objective
has several local optima: a break-point term, a nonlinear term, a
score-driven filter or a latent-state mixture.

## Usage

``` r
start_search(optimizer = optimizers7::sa(), over = NULL)
```

## Arguments

- optimizer:

  The optimizer to search with, any optimizers7 optimizer. Defaults to
  [`optimizers7::sa()`](https://statmodels7.github.io/optimizers7/reference/sa.html),
  simulated annealing, which is the one built for a multimodal surface.

- over:

  Optional character vector naming the coefficients to search over,
  overriding the choice above. Names not in the model are an error.

## Value

A `StartSearch` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md),
with properties `optimizer` and `over`.

## Once, before the fit

A search belongs to the starting value, never to the scoring step. Given
to `inner_optimizer` instead, as `chain(sa(), iwls())`, it would rerun
at every hyperparameter the outer criterion tried, which on an ordinary
fit is 46 times and inside
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md) is
the folds times the path, returning the same answer each time.

## On the likelihood alone

The penalties are off. What a starting value has to get right is the
basin of the likelihood, the one thing the fit will not correct by
itself. The penalties enter afterwards, when their hyperparameters are
estimated, and at the probe values they represent nobody's choice: a
lasso at \\\lambda = 1\\ against an unaveraged log-likelihood empties
whole blocks, so a search there would explore a different model.

## Which coordinates are searched

The non-convex ones, which the toolkit can already name: the parameters
of a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)),
the coefficients of a term that recomputes its own block
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
and the break-point terms), and each equation's intercept.

Everything else keeps its default. A smooth, a ridge or a random effect
is a convex block whose optimum the scoring step reaches from anywhere,
and searching over a thousand random-effect coefficients would spend the
whole budget on the part of the model that needs it least.

A penalized coordinate is excluded **wherever it sits**, including
inside a non-convex term. A sub-formula develops a break-point or a
nonlinear parameter over groups, and those deviations are columns of the
term's own block. They are the case the rule exists for: on the
likelihood alone nothing identifies them, since the penalty is what
does, so a search over them fits each group's own points and moves away
from the penalized mode. Measured on
`jseg(x, npsi = 2, by = ~random(~1|id))` over thirty groups, 210 of the
219 coordinates are such deviations.

`over` overrides the choice by name.

## The hyperparameters are not searched

They cannot be, from here. The objective is the likelihood with the
penalties off, in which a hyperparameter does not appear, so no proposal
could change one.

A global search over them is a search over the **outer** criterion,
where each proposal costs a full refit instead of one likelihood
evaluation, and that is `statmod(outer_optimizer = optimizers7::sa())`.
A kinked penalty is outside even that: its hyperparameter has a known
upper end and is swept by a warm-started path, which a random jump would
fail to improve on and would destroy.

## See also

[`start_intercepts()`](https://statmodels7.github.io/statmodels7/reference/start_intercepts.md)
for the default,
[`optimizers7::sa()`](https://statmodels7.github.io/optimizers7/reference/sa.html)
for the search,
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
for the `outer_optimizer` argument that searches the hyperparameters
instead.

## Examples

``` r
start_search()
#> <start> search with simulated annealing (uniform)
start_search(optimizers7::sa(maxit = 20))
#> <start> search with simulated annealing (uniform)

# It changes nothing on a convex fit: there are no non-convex coordinates
# to search but the intercepts, and the scoring step reaches those anyway.
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.3)
a <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
b <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd,
             start = start_search(optimizers7::sa(maxit = 30)))
max(abs(unlist(a@coefficients) - unlist(b@coefficients)))
#> [1] 0
```
