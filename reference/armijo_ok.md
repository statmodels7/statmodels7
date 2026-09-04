# Armijo's Sufficient-Decrease Condition

Whether a trial point may be accepted: the objective has fallen by at
least a fixed share of what the model predicted over that step.

## Usage

``` r
armijo_ok(vnew, value, step, gd, c1 = 1e-04)
```

## Arguments

- vnew:

  The objective at the trial point.

- value:

  The objective at the current one.

- step:

  The step length.

- gd:

  \\g^\top\delta\\, negative for a descent direction.

- c1:

  The share of the predicted decrease that is asked for.

## Value

A single logical.

## Details

\$\$f(x + s d) \\\le\\ f(x) + c_1 s\\ g^\top d .\$\$ **The sign of
\\g^\top d\\ is the whole of it.** The increment solves \\(H + S)\delta
= -g\\, so \\g^\top\delta = -g^\top(H+S)^{-1}g\\ is NEGATIVE – measured,
in 19 of 19 solves of one fit, over \\\[-6435, -0.0119\]\\ – and the
bound therefore sits BELOW \\f(x)\\ and asks for a decrease. Written
with the term SUBTRACTED it changes sign and the test permits an
INCREASE of \\c_1 s \|g^\top\delta\|\\, which on that fit was as much as
0.6435 a step, where the loop's own comment said sufficient decrease
rather than mere non-increase. It is named rather than written inline so
that the sign carries a test of its own: a test over a fit sees the
condition only through the steps that fit happens to take, and none of
them relied on the extra licence.

## See also

[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md),
the loop that asks it.
