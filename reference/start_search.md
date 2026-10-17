# Search the Likelihood for a Starting Point

Runs a global search once, before the fit begins, over the coordinates
where the likelihood is not convex. The strategy is built for a model
whose objective has several local optima: a break-point term, a
nonlinear term, a score-driven filter or a latent-state mixture.

## Usage

``` r
start_search(optimizer = optimizers7::sa(), over = NULL)
```

## Arguments

- optimizer:

  The optimizer that performs the search, any optimizers7 optimizer.
  Defaults to
  [`optimizers7::sa()`](https://statmodels7.github.io/optimizers7/reference/sa.html),
  simulated annealing, which is the one built for a multimodal surface.

- over:

  Optional character vector naming the coefficients to search over,
  overriding the choice above. A coefficient is named by its label, the
  parameter and the coefficient joined by a colon
  (`"sigma:(Intercept)"`); a bare name is accepted where only one
  equation carries it. Names not in the model, and a bare name that
  several equations carry, are an error.

## Value

A `StartSearch` object, inheriting from
[`start_strategy()`](https://statmodels7.github.io/statmodels7/reference/start_strategy.md),
with properties `optimizer` and `over`.

## Once, before the fit

A search is a starting value and is not part of the scoring step. If it
were given to `inner_optimizer`, it would be rerun at every
hyperparameter value that the outer criterion tries and, inside
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md), at
every point of the path for every fold, and each run would reach the
same point.

## On the likelihood alone

The penalties are off. A starting value has to select the basin of the
likelihood, which is the one thing that the fit does not correct by
itself. The penalties enter afterwards, when their hyperparameters are
estimated, and at the probe values of the hyperparameters they are not a
chosen penalty, so a search with them would explore the likelihood of a
different model.

## Which coordinates are searched

The search covers the non-convex coordinates, which the toolkit can
name: the parameters of a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)),
the coefficients of a term that recomputes its own block
([`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html),
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html)
and the break-point terms), and each equation's intercept.

All other coefficients keep their default values. A smooth, a ridge or a
random effect is a convex block whose optimum the scoring step reaches
from any start, and a search over a thousand random-effect coefficients
would spend the budget on the part of the model that needs a search
least.

A penalized coordinate is excluded wherever it sits, including inside a
non-convex term. A sub-formula can develop a break-point or a nonlinear
parameter over groups, and the group deviations are columns of the block
of that term. Without the penalty these deviations are not identified,
because only the penalty identifies them, so a search over them would
fit the points of each group and move away from the penalized mode. In a
`jseg()` term with `by = ~random(~1|id)`, most coordinates are such
deviations.

The argument `over` replaces this choice by naming the coefficients. A
structural term's own parameters are searched in either case.

## The hyperparameters

The search does not cover the hyperparameters. Its objective is the
likelihood with the penalties off, in which the hyperparameters do not
appear, so a proposal cannot change them.

A global search over the hyperparameters is a search over the outer
criterion, where each proposal requires a full refit and not one
likelihood evaluation. It is requested by
`statmod(outer_optimizer = optimizers7::sa())`. The hyperparameters of a
kinked penalty are outside that search as well: they are chosen along a
warm-started path scored by `sparse_criterion`, which a random jump
would not improve.

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
