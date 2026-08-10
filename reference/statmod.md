# Fit a Model

Reads one formula carrying every parameter of a distribution, assembles
the terms it names into a penalized likelihood, and fits it.

## Usage

``` r
statmod(
  formula,
  distrib,
  data,
  weights = NULL,
  offsets = NULL,
  inner_method = iwls(),
  start = NULL,
  maxit = 50L,
  tol = 1e-08,
  verbose = 0
)
```

## Arguments

- formula:

  The model formula.

- distrib:

  A distributions7 distribution object.

- data:

  A data frame.

- weights:

  Optional prior weights, taken as given and not normalized.

- offsets:

  Optional named list of offsets, one per parameter.

- inner_method:

  How the smooth block is fitted:
  [`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
  or an optimizers7 optimizer.

- start:

  Optional starting coefficients, a named list.

- maxit:

  The alternation's iteration budget.

- tol:

  The alternation's tolerance, on the relative change in the objective.

- verbose:

  A level from 0 to 3, or a named logical vector.

## Value

An object of class
[`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Details

**The formula.** The equations of the distribution's parameters are
separated by `|`, the first carrying the response:

        y ~ x1 + ridge(R) + lasso(L)  |  sigma ~ z  |  nu ~ 1

A parameter with no equation gets an intercept. See
[`statmod_equations`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md),
whose recovery is not the obvious one.

**The fitting scheme.** The terms split in two by a property each one
already reports. Every term whose penalty is twice differentiable in its
coefficients – an unpenalized block, a ridge, a spline, a random effect
– is estimated in ONE system by `inner_method`, because their joint
curvature exists and using it is what makes a fit converge in a handful
of iterations. A term whose penalty has a kink – lasso, scad, mcp – is
estimated by a method of its own with everything else held fixed. The
fit alternates between the two until the objective and every block stop
moving.

**The objective is unaveraged**: minus the weighted log-likelihood plus
the penalties at full size, since a penalty is a negative log-prior and
a posterior adds the two at full size. What is scaled instead is the
stopping rule, so that a threshold means the same thing at \\n = 10\\
and at \\n = 10^7\\.

**Verbosity** has three levels, naming the loops rather than counting
them: `1` the alternation, `2` the inner method's own iterations, `3`
the optimizers' traces as well. A named form is accepted too, as
`verbose = c(blocks = TRUE, inner = TRUE, optimizer = FALSE)`, since
watching the alternation while silencing a chatty inner optimizer is the
common case.

## See also

[`statmod_spec`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
[`iwls`](https://statmodels7.github.io/statmodels7/reference/iwls.md),
[`loglik`](https://statmodels7.github.io/statmodels7/reference/loglik.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(60))
dd$y <- 1 + 2 * dd$x + rnorm(60, sd = 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
fit
#> <statmodels7::StatmodFit>
#>  @ spec        : <statmodels7::StatmodSpec>
#>  .. @ formula   :Class 'formula'  language y ~ x
#>  .. ..  ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. @ distrib   : <distributions7::Gaussian1Distrib>
#>  .. .. @ distrib_name         : chr "gaussian1"
#>  .. .. @ dimension            : chr "univariate"
#>  .. .. @ bounds               : num [1:2] -Inf Inf
#>  .. .. @ params               : chr [1:2] "mu" "sigma"
#>  .. .. @ params_interpretation: Named chr [1:2] "mean" "standard deviation"
#>  .. .. .. - attr(*, "names")= chr [1:2] "mu" "sigma"
#>  .. .. @ n_params             : num 2
#>  .. .. @ params_bounds        :List of 2
#>  .. .. .. $ mu   : num [1:2] -Inf Inf
#>  .. .. .. $ sigma: num [1:2] 0 Inf
#>  .. .. @ link_params          :List of 2
#>  .. .. .. $ mu   : <linkfunctions7::IdentityLink>
#>  .. .. ..  ..@ link_name  : chr "identity"
#>  .. .. ..  ..@ link_bounds: num [1:2] -Inf Inf
#>  .. .. ..  ..@ link_params: NULL
#>  .. .. .. $ sigma: <linkfunctions7::LogLink>
#>  .. .. ..  ..@ link_name  : chr "log"
#>  .. .. ..  ..@ link_bounds: num [1:2] 0 Inf
#>  .. .. ..  ..@ link_params: NULL
#>  .. .. @ params_smooth        : logi(0) 
#>  .. @ equations :List of 2
#>  .. .. $ mu   :Class 'formula'  language ~x
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. .. $ sigma:Class 'formula'  language ~1
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. @ terms     :List of 2
#>  .. .. $ mu   :List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:60, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:60] "1" "2" "3" "4" ...
#>  .. .. .. .. ..  ..$ : chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ coef_names: chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ blueprint :List of 3
#>  .. .. .. .. .. $ terms    :Classes 'terms', 'formula'  language ~x
#>  .. .. .. .. ..  .. ..- attr(*, "variables")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "factors")= int [1, 1] 1
#>  .. .. .. .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "term.labels")= chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "order")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "response")= int 0
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e004b98c00> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr "numeric"
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr "x"
#>  .. .. .. .. .. $ xlev     : Named list()
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~x
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e004b98c00> 
#>  .. .. $ sigma:List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:60, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:60] "1" "2" "3" "4" ...
#>  .. .. .. .. ..  ..$ : chr "(Intercept)"
#>  .. ..  .. ..@ coef_names: chr "(Intercept)"
#>  .. ..  .. ..@ blueprint :List of 3
#>  .. .. .. .. .. $ terms    :Classes 'terms', 'formula'  language ~1
#>  .. .. .. .. ..  .. ..- attr(*, "variables")= language list()
#>  .. .. .. .. ..  .. ..- attr(*, "factors")= int(0) 
#>  .. .. .. .. ..  .. ..- attr(*, "term.labels")= chr(0) 
#>  .. .. .. .. ..  .. ..- attr(*, "order")= int(0) 
#>  .. .. .. .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "response")= int 0
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e004b98c00> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list()
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr(0) 
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr(0) 
#>  .. .. .. .. .. $ xlev     : NULL
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~1
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e004b98c00> 
#>  .. @ response  : num [1:60] 2.21 1.693 2.34 2.79 0.715 ...
#>  .. @ n_obs     : int 60
#>  .. @ weights   : num [1:60] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. @ offsets   :List of 2
#>  .. .. $ mu   : NULL
#>  .. .. $ sigma: NULL
#>  .. @ intercepts: Named logi [1:2] TRUE TRUE
#>  .. .. - attr(*, "names")= chr [1:2] "mu" "sigma"
#>  @ coefficients:List of 2
#>  .. $ mu   : num [1:2] 1.09 1.95
#>  .. $ sigma: num -0.836
#>  @ hyper       :List of 2
#>  .. $ mu   : list()
#>  .. $ sigma: list()
#>  @ loglik      : num -34.9
#>  @ objective   : num 34.9
#>  @ edf         :'data.frame':    2 obs. of  4 variables:
#>  .. $ parameter   : chr  "mu" "sigma"
#>  .. $ term        : chr  "linpar" "linpar"
#>  .. $ coefficients: int  2 1
#>  .. $ edf         : num  2 1
#>  @ fitted      :List of 2
#>  .. $ mu   : num [1:60] 1.61 1.81 2.2 2.86 1.48 ...
#>  .. $ sigma: num [1:60] 0.433 0.433 0.433 0.433 0.433 ...
#>  @ converged   : logi TRUE
#>  @ elapsed     : num 0.026
#>  @ history     :List of 3
#>  .. $ outer : NULL
#>  .. $ blocks:'data.frame':   1 obs. of  6 variables:
#>  ..  ..$ sweep     : int 1
#>  ..  ..$ block     : chr "smooth"
#>  ..  ..$ objective : num 34.9
#>  ..  ..$ change    : num 164
#>  ..  ..$ iterations: int 11
#>  ..  ..$ converged : logi TRUE
#>  .. $ inner :'data.frame':   10 obs. of  8 variables:
#>  ..  ..$ iteration: int [1:10] 1 2 3 4 5 6 7 8 9 10
#>  ..  ..$ objective: num [1:10] 169.4 139.8 110.7 83.1 59.1 ...
#>  ..  ..$ score    : num [1:10] 3.805 0.996 0.989 0.97 0.92 ...
#>  ..  ..$ step     : num [1:10] 1 1 1 1 1 1 1 1 1 1
#>  ..  ..$ rank     : int [1:10] 3 3 3 3 3 3 3 3 3 3
#>  ..  ..$ route    : chr [1:10] "qr" "qr" "qr" "qr" ...
#>  ..  ..$ sweep    : int [1:10] 1 1 1 1 1 1 1 1 1 1
#>  ..  ..$ block    : chr [1:10] "smooth" "smooth" "smooth" "smooth" ...
#>  @ methods     :List of 2
#>  .. $ smooth: <statmodels7::Iwls>
#>  ..  ..@ hessian      : chr "expected"
#>  ..  ..@ approx       : chr "bartlett"
#>  ..  ..@ decomposition: chr "qr"
#>  ..  ..@ maxit        : num 100
#>  ..  ..@ tol          : num 1e-06
#>  ..  ..@ step_halving : num 30
#>  .. $ sparse: chr(0) 
#>  @ call        : language statmod(formula = y ~ x, distrib = distributions7::gaussian1_distrib(),      data = dd)

# every parameter can be modelled
statmod(y ~ x | sigma ~ x, distributions7::gaussian1_distrib(), dd)
#> <statmodels7::StatmodFit>
#>  @ spec        : <statmodels7::StatmodSpec>
#>  .. @ formula   :Class 'formula'  language y ~ x | sigma ~ x
#>  .. ..  ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. @ distrib   : <distributions7::Gaussian1Distrib>
#>  .. .. @ distrib_name         : chr "gaussian1"
#>  .. .. @ dimension            : chr "univariate"
#>  .. .. @ bounds               : num [1:2] -Inf Inf
#>  .. .. @ params               : chr [1:2] "mu" "sigma"
#>  .. .. @ params_interpretation: Named chr [1:2] "mean" "standard deviation"
#>  .. .. .. - attr(*, "names")= chr [1:2] "mu" "sigma"
#>  .. .. @ n_params             : num 2
#>  .. .. @ params_bounds        :List of 2
#>  .. .. .. $ mu   : num [1:2] -Inf Inf
#>  .. .. .. $ sigma: num [1:2] 0 Inf
#>  .. .. @ link_params          :List of 2
#>  .. .. .. $ mu   : <linkfunctions7::IdentityLink>
#>  .. .. ..  ..@ link_name  : chr "identity"
#>  .. .. ..  ..@ link_bounds: num [1:2] -Inf Inf
#>  .. .. ..  ..@ link_params: NULL
#>  .. .. .. $ sigma: <linkfunctions7::LogLink>
#>  .. .. ..  ..@ link_name  : chr "log"
#>  .. .. ..  ..@ link_bounds: num [1:2] 0 Inf
#>  .. .. ..  ..@ link_params: NULL
#>  .. .. @ params_smooth        : logi(0) 
#>  .. @ equations :List of 2
#>  .. .. $ mu   :Class 'formula'  language ~x
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. .. $ sigma:Class 'formula'  language ~x
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e0031680e8> 
#>  .. @ terms     :List of 2
#>  .. .. $ mu   :List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:60, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:60] "1" "2" "3" "4" ...
#>  .. .. .. .. ..  ..$ : chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ coef_names: chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ blueprint :List of 3
#>  .. .. .. .. .. $ terms    :Classes 'terms', 'formula'  language ~x
#>  .. .. .. .. ..  .. ..- attr(*, "variables")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "factors")= int [1, 1] 1
#>  .. .. .. .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "term.labels")= chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "order")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "response")= int 0
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e00577ab80> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr "numeric"
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr "x"
#>  .. .. .. .. .. $ xlev     : Named list()
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~x
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e00577ab80> 
#>  .. .. $ sigma:List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:60, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:60] "1" "2" "3" "4" ...
#>  .. .. .. .. ..  ..$ : chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ coef_names: chr [1:2] "(Intercept)" "x"
#>  .. ..  .. ..@ blueprint :List of 3
#>  .. .. .. .. .. $ terms    :Classes 'terms', 'formula'  language ~x
#>  .. .. .. .. ..  .. ..- attr(*, "variables")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "factors")= int [1, 1] 1
#>  .. .. .. .. ..  .. .. ..- attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. .. .. ..$ : chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "term.labels")= chr "x"
#>  .. .. .. .. ..  .. ..- attr(*, "order")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "intercept")= int 1
#>  .. .. .. .. ..  .. ..- attr(*, "response")= int 0
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e00577ab80> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr "numeric"
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr "x"
#>  .. .. .. .. .. $ xlev     : Named list()
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~x
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e00577ab80> 
#>  .. @ response  : num [1:60] 2.21 1.693 2.34 2.79 0.715 ...
#>  .. @ n_obs     : int 60
#>  .. @ weights   : num [1:60] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. @ offsets   :List of 2
#>  .. .. $ mu   : NULL
#>  .. .. $ sigma: NULL
#>  .. @ intercepts: Named logi [1:2] TRUE TRUE
#>  .. .. - attr(*, "names")= chr [1:2] "mu" "sigma"
#>  @ coefficients:List of 2
#>  .. $ mu   : num [1:2] 1.11 1.91
#>  .. $ sigma: num [1:2] -0.508 -0.691
#>  @ hyper       :List of 2
#>  .. $ mu   : list()
#>  .. $ sigma: list()
#>  @ loglik      : num -33.4
#>  @ objective   : num 33.4
#>  @ edf         :'data.frame':    2 obs. of  4 variables:
#>  .. $ parameter   : chr  "mu" "sigma"
#>  .. $ term        : chr  "linpar" "linpar"
#>  .. $ coefficients: int  2 2
#>  .. $ edf         : num  2 2
#>  @ fitted      :List of 2
#>  .. $ mu   : num [1:60] 1.62 1.82 2.2 2.84 1.49 ...
#>  .. $ sigma: num [1:60] 0.501 0.465 0.405 0.321 0.524 ...
#>  @ converged   : logi TRUE
#>  @ elapsed     : num 0.046
#>  @ history     :List of 3
#>  .. $ outer : NULL
#>  .. $ blocks:'data.frame':   1 obs. of  6 variables:
#>  ..  ..$ sweep     : int 1
#>  ..  ..$ block     : chr "smooth"
#>  ..  ..$ objective : num 33.4
#>  ..  ..$ change    : num 166
#>  ..  ..$ iterations: int 19
#>  ..  ..$ converged : logi TRUE
#>  .. $ inner :'data.frame':   18 obs. of  8 variables:
#>  ..  ..$ iteration: int [1:18] 1 2 3 4 5 6 7 8 9 10 ...
#>  ..  ..$ objective: num [1:18] 169.9 141.6 115.5 93.4 75.4 ...
#>  ..  ..$ score    : num [1:18] 3.805 0.978 0.945 0.87 0.737 ...
#>  ..  ..$ step     : num [1:18] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  ..$ rank     : int [1:18] 4 4 4 4 4 4 4 4 4 4 ...
#>  ..  ..$ route    : chr [1:18] "qr" "qr" "qr" "qr" ...
#>  ..  ..$ sweep    : int [1:18] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  ..$ block    : chr [1:18] "smooth" "smooth" "smooth" "smooth" ...
#>  @ methods     :List of 2
#>  .. $ smooth: <statmodels7::Iwls>
#>  ..  ..@ hessian      : chr "expected"
#>  ..  ..@ approx       : chr "bartlett"
#>  ..  ..@ decomposition: chr "qr"
#>  ..  ..@ maxit        : num 100
#>  ..  ..@ tol          : num 1e-06
#>  ..  ..@ step_halving : num 30
#>  .. $ sparse: chr(0) 
#>  @ call        : language statmod(formula = y ~ x | sigma ~ x, distrib = distributions7::gaussian1_distrib(),      data = dd)
```
