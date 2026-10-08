# The Interval and Standard Deviation of a New Group's Parameter

For `interval = "group"`: the distribution of a new group's parameter
\\\theta^\* = h^{-1}(\eta^\*)\\, with \\\eta^\* = \eta_0 + e + z^\top
b\\, where \\e\\ is the estimation error of the fixed part and \\b\\ the
effects of a new group, drawn from the prior the fit estimated. Returns
the ends of its interval and its standard deviation, on the predictor's
scale and on the parameter's.

## Usage

``` r
group_interval(object, spec, design, aside, ep, su, level, ...)
```

## Arguments

- object:

  The
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- spec, design:

  The specification and design at the rows predicted.

- aside:

  The rows of
  [`random_modes()`](https://statmodels7.github.io/statmodels7/reference/random_modes.md)
  that are not `"conditional"`.

- ep:

  The predictors at the fixed part, as
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md)
  returns them.

- su:

  The standard errors of
  [`predict_se()`](https://statmodels7.github.io/statmodels7/reference/predict_se.md).

- level:

  The interval's level.

- ...:

  Passed to
  [`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md).

## Value

`su`, with `se_eta`, `eta_lower`, `eta_upper`, `se`, `lower` and `upper`
replaced for every parameter a new group's effects reach.

## Details

The inverse link is monotone, so a quantile of \\\theta^\*\\ is the
inverse link of the same quantile of \\\eta^\*\\, and the interval is
exact whatever the link. Under Gaussian priors \\\eta^\*\\ is Gaussian
with variance \\s^2 = \mathrm{se}\_0^2 + z^\top\Sigma_b z\\, and the
ends are \\h^{-1}(\eta_0 \pm z\_{(1+\ell)/2}\\ s)\\. Under any other
prior \\\eta^\*\\ is a mixture of Gaussians: a scale mixture over the
nodes of
[`gamma_nodes()`](https://statmodels7.github.io/statmodels7/reference/gamma_nodes.md)
for a Student t, the prior's own quantiles for one univariate prior of
another family
([`group_quantile_nodes()`](https://statmodels7.github.io/statmodels7/reference/group_quantile_nodes.md)),
and the Monte Carlo draws of
[`predictive_mixture()`](https://statmodels7.github.io/statmodels7/reference/predictive_mixture.md)
for anything else. Each end is found by bisection on the mixture's
distribution function.

The standard deviation is that of \\\theta^\*\\, not the delta method,
which under a log link and \\s = 0.85\\ is 42 per cent too small. Under
a Gaussian \\\eta^\*\\ it is \\s\\ for the identity link, \\e^{\eta_0 +
s^2/2}\sqrt{e^{s^2} - 1}\\ for the log link and a 40-node Gauss-Hermite
rule for any other link. It is reported only where it is known to exist:

- never where the predictor's domain is not the whole line (the square
  root and inverse links), since a Gaussian \\\eta^\*\\ leaves it;

- always for an inverse link bounded on both sides (logit, probit,
  cloglog, a bounded link), from the mixture's nodes;

- for Gaussian priors under every other link;

- under a prior that is not Gaussian and an unbounded link, only for the
  identity link, where the variance is \\\mathrm{se}\_0^2 +
  z^\top\mathrm{Var}(b)\\z\\: a Student t gives \\\Sigma\nu/(\nu - 2)\\,
  `NA` for \\\nu \le 2\\, and a univariate prior of another family its
  own variance. Under a log link a Student t has no moment generating
  function, so \\E\[\theta^\*\]\\ is infinite, and the standard
  deviation is `NA`.
