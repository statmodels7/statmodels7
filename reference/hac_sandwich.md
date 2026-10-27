# HAC Sandwich Variance of the Coefficients

Select a heteroskedasticity and autocorrelation consistent meat for the
frequentist variance \\V = V_b B V_b\\, with \\V_b = (H + S)^{-1}\\ the
inverse of the penalized information. Pass the object as the `sandwich`
argument of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).
It is the counterpart of
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
for dependent observations, and the estimator of Newey and West (1987)
with the Bartlett kernel.

## Usage

``` r
hac_sandwich(
  kernel = c("bartlett", "parzen", "truncated"),
  bandwidth = NULL,
  time = NULL,
  cluster = NULL,
  data = NULL
)
```

## Arguments

- kernel:

  `"bartlett"` (the default), `"parzen"` or `"truncated"`.

- bandwidth:

  `NULL` (the default), which takes \\\lfloor 4 (T/100)^{2/9} \rfloor\\,
  or a non-negative whole number, the largest lag \\L\\.

- time:

  `NULL` (the default) or a one-sided formula with one variable, such as
  `~ t`, for a model without a structural term.

- cluster:

  `NULL` (the default) or a one-sided formula with one variable, such as
  `~ id`, that separates the series of a model without a structural
  term.

- data:

  `NULL` (the default) or the data frame that holds `time` and
  `cluster`, with one row for each observation of the fit.

## Value

A `HacSandwich` object with the properties `type` (`"HAC"`), `kernel`,
`bandwidth`, `time`, `cluster` and `data`.

## Details

The score \\u_t\\ of observation \\t\\ is the derivative of its
log-likelihood contribution in the whole vector of unknowns: the
coefficients of every distribution parameter, and the own parameters of
a structural term on the unconstrained scale. For a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term the contribution of \\t\\ depends on the parameters through the
filtered level, and \\u_t\\ carries that dependence: it is the sum over
the equations of the derivative of the log-density in the predictor, on
the link scale, times the derivative of the predictor in the unknowns,
which the forward recursion of the filter supplies. For a
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term the contribution is the increment \\\log p(y_t \mid y_1, \dots,
y\_{t-1})\\ of the forward recursion, and \\u_t\\ its derivative,
computed in forward mode through the recursion. In both cases \\\sum_t
u_t\\ is the gradient of the log-likelihood.

With a well specified filter the scores are martingale differences and
\\B = \sum_t u_t u_t'\\ already estimates their variance. The HAC meat
adds the autocovariances, \$\$B = \Gamma_0 + \sum\_{l=1}^{L}
k\\\left(\frac{l}{L + 1}\right) (\Gamma_l + \Gamma_l'), \qquad \Gamma_l
= \sum_t u_t u\_{t-l}',\$\$ where \\k\\ is the kernel and \\L\\ the
bandwidth, the largest lag. The Bartlett kernel is \\k(x) = 1 - x\\, the
Parzen kernel is \\1 - 6x^2 + 6x^3\\ for \\x \le 1/2\\ and \\2(1 -
x)^3\\ up to one, and the truncated kernel is one. With `bandwidth = 0`
the meat is \\\sum_t u_t u_t'\\, the meat of
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md).
No finite-sample correction is applied.

The default bandwidth is \\L = \lfloor 4 (T / 100)^{2/9} \rfloor\\, with
\\T\\ the number of observations of a group (the mean number, when the
groups differ in length). A bandwidth larger than the longest group
minus one is reduced to that.

With a `by` in the structural term the series are panels. A lag runs
within a group only, so the product of two scores of different groups
enters no \\\Gamma_l\\, and the groups are independent. Observations are
ordered within a group as the structural term orders them (its `time`),
and a lag is a step in that order, whatever the gaps between times.

A model without a structural term has no order of its own. `time` names
the variable that orders the observations and `cluster` the variable
that separates independent series; both are read from `data`, or from
the data of the call to
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
which is evaluated again in the environment of the model formula. For a
model with a structural term they are not accepted, since the term fixes
the order and the groups.

The prior weights of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
multiply the log-likelihood contribution of an observation, and so its
score; they are not otherwise interpreted. Models with a term whose
block moves with its coefficients (the break-point terms and
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html))
are rejected.

## References

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
heteroskedasticity and autocorrelation consistent covariance matrix.
*Econometrica*, 55(3), 703–708.

## See also

[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md),
which takes the object, and
[`hc0_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hc0_sandwich.md)
for independent observations.

## Examples

``` r
set.seed(1)
n <- 120
dd <- data.frame(t = seq_len(n), x = rnorm(n))
e <- as.numeric(stats::filter(rnorm(n), 0.6, method = "recursive"))
dd$y <- 1 + 0.5 * dd$x + e
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hc0_sandwich()))))[1:2]
#> mu:(Intercept)           mu:x 
#>      0.1152153      0.1236192 
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hac_sandwich(time = ~ t)))))[1:2]
#> mu:(Intercept)           mu:x 
#>      0.1894405      0.1017458 
```
