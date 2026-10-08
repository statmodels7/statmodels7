# Test That a Break-Point Exists

Tests the hypothesis that a model has no break-point against the
alternative of the break-point term it carries: no change of slope for a
[`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
no change of level for a
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html),
and neither for a
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html).
The term can sit in the equation of any distribution parameter, beside
any other term, and inside a parameter of
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html).

## Usage

``` r
statmod_breakpoint_test(
  fit,
  param = NULL,
  term = NULL,
  k = 10L,
  n_boot = 199L,
  seed = NULL
)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md)
  carrying the break-point term.

- param:

  The distribution parameter whose equation carries the term, or `NULL`
  to find it.

- term:

  The label of the term (`"seg"`, `"jump"`, `"jseg"` by default), or
  `NULL` to find it. Needed only where the model carries more than one
  break-point term.

- k:

  The number of positions at which a continuous process is read.

- n_boot:

  The number of bootstrap replicates for a sharp step.

- seed:

  A seed for the bootstrap, or `NULL`.

## Value

An object of class `"htest"` with the supremum of the score process as
`statistic`, the number of change coefficients as `parameter`, the
`p.value`, the position of the supremum as `estimate`, and the process
itself as `process`, a data frame with columns `psi` and `score`.

## Details

Under the null hypothesis the position \\\psi\\ of the break-point does
not enter the model, so the classical tests of the change against zero
have no fixed reference distribution: the position is a nuisance
parameter present only under the alternative (Davies, 1987). The test
reads Rao's score statistic for the change at a fixed position \\q\\,
\$\$S(q) = U_c(q)' \\ \[K(q)^{-1}\]\_{cc} \\ U_c(q),\$\$ where
\\U_c(q)\\ is the score of the change coefficients at the fit of the
null model and \\K(q)\\ the penalized information of the model with the
break-point held at \\q\\. One fit of the null model serves every \\q\\:
the process costs one gradient and one Hessian per position. The
statistic of the test is \\M = \sup_q S(q)\\.

Where the position moves the model continuously (a change of slope, or a
smoothed step), \\S(q)\\ is evaluated at `k` positions equally spaced
between the confinement limits of the term, and the p-value is the upper
bound of Davies (1987) for a \\\chi^2_s\\ process, \$\$P(\sup_q S(q) \>
M) \le P(\chi^2_s \> M) + V M^{(s-1)/2} e^{-M/2} 2^{-s/2} /
\Gamma(s/2),\$\$ with \\V = \sum_j \|S(q\_{j+1})^{1/2} -
S(q_j)^{1/2}\|\\ the total variation of the root of the process and
\\s\\ the number of change coefficients. This is the test of
[`segmented::davies.test()`](https://rdrr.io/pkg/segmented/man/davies.test.html),
read on the score of the model at hand.

Where the step is sharp the process is constant between consecutive
values of the covariate and jumps at each of them, and the bound, which
counts the variation over every interval, is too wide to be useful. The
supremum is then taken over every interval inside the confinement
limits, and its distribution under the null is estimated by a parametric
bootstrap: `n_boot` responses are drawn from the null fit, the null
model is refitted to each, and the p-value is the share of replicates
whose supremum reaches the observed one, \\(1 + \\\\M^\* \ge M\\)/(1 +
B)\\. The draws are conditional on the fitted random effects, if any.

The hyperparameters are held at the values of the fit, so the test is
conditional on the smoothing the data chose. The term must carry one
break-point whose position has no development; for a choice between one
and several break-points, compare the fits by an information criterion.

## References

Davies, R. B. (1987). Hypothesis testing when a nuisance parameter is
present only under the alternative. *Biometrika*, 74, 33–43.

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22, 3055–3071.

## See also

[`statmod_test()`](https://statmodels7.github.io/statmodels7/reference/statmod_test.md)
for a test of one coefficient against one value.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(120, 0, 10))
dd$y <- 1 + 0.5 * dd$x - 0.8 * pmax(dd$x - 6, 0) + rnorm(120, sd = 0.5)
fit <- statmod(y ~ seg(x), distributions7::gaussian1_distrib(), dd)
statmod_breakpoint_test(fit)
#> 
#>  Davies' test for a break-point (score process, upper bound)
#> 
#> data:  seg in the equation of mu
#> sup S = 69.06, df = 1, p-value = 3.549e-15
#> alternative hypothesis: a change of slope
#> sample estimates:
#>      psi 
#> 6.187742 
#> 
```
