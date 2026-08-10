# Simulate a Response From a Written Model

Takes a formula, a distribution and a data frame of covariates, and
draws a response from that model – at coefficients the caller supplies,
or at coefficients drawn at random.

## Usage

``` r
rstatmod(formula, distrib, data, par = NULL, sd = 1, offsets = NULL)
```

## Arguments

- formula:

  The model formula, as
  [`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
  takes it.

- distrib:

  A distributions7 distribution object.

- data:

  A data frame of covariates.

- par:

  Optional named list of coefficient vectors.

- sd:

  The standard deviation of the drawn coefficients.

- offsets:

  Optional named list of offsets.

## Value

The data frame with the response added, carrying the attributes `"par"`
(the coefficients used) and `"theta"` (the parameters they gave).

## Details

This is not [`simulate`](https://rdrr.io/r/stats/simulate.html), which
draws from a model that has already been fitted. The `r` prefix is R's
own for a random draw, so the two names cannot be confused.

The point of it is to have data whose truth is known: write the model,
draw from it, fit it back, and see whether the fit recovers what was put
in. A covariate needs no declaring – a factor becomes its contrasts and
a numeric stays itself, because the design comes from the same
interpreter a fit uses.

**The coefficients.** `par = NULL` draws them, each independently from
`rnorm(1, 0, sd)`, which on the link scale gives predictors of order
one. A named list fixes them instead, one vector per distribution
parameter in the design's order; a parameter left out of that list is
drawn. [`coef()`](https://rdrr.io/r/stats/coef.html) on the result
reports what was used, drawn or given, so a simulation is reproducible
from its own output.

**The response's name** is the formula's left-hand side when it is a
symbol, and `"y"` otherwise.

## See also

[`statmod`](https://statmodels7.github.io/statmodels7/reference/statmod.md),
[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(50), g = factor(rep(c("a", "b"), 25)))

# coefficients drawn
sim <- rstatmod(y ~ x + g, distributions7::gaussian1_distrib(), dd)
attr(sim, "par")
#> $mu
#> (Intercept)           x          gb 
#> -0.05612874 -0.15579551 -1.47075238 
#> 
#> $sigma
#> (Intercept) 
#>  -0.4781501 
#> 

# or given, and recovered by a fit
sim2 <- rstatmod(y ~ x, distributions7::gaussian1_distrib(), dd,
                 par = list(mu = c(1, 2), sigma = log(0.3)))
statmod(y ~ x, distributions7::gaussian1_distrib(), sim2)
#> <statmodels7::StatmodFit>
#>  @ spec        : <statmodels7::StatmodSpec>
#>  .. @ formula   :Class 'formula'  language y ~ x
#>  .. ..  ..- attr(*, ".Environment")=<environment: 0x55dffe86bcb0> 
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
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55dffe86bcb0> 
#>  .. .. $ sigma:Class 'formula'  language ~1
#>  .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55dffe86bcb0> 
#>  .. @ terms     :List of 2
#>  .. .. $ mu   :List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:50, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:50] "1" "2" "3" "4" ...
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
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e00116ca68> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list(x)
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr "numeric"
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr "x"
#>  .. .. .. .. .. $ xlev     : Named list()
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~x
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e00116ca68> 
#>  .. .. $ sigma:List of 1
#>  .. ..  ..$ linpar: <modelterms7::LinparTerm>
#>  .. ..  .. ..@ label     : chr ""
#>  .. ..  .. ..@ X         : num [1:50, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. .. .. .. .. - attr(*, "dimnames")=List of 2
#>  .. .. .. .. ..  ..$ : chr [1:50] "1" "2" "3" "4" ...
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
#>  .. .. .. .. ..  .. ..- attr(*, ".Environment")=<environment: 0x55e00116ca68> 
#>  .. .. .. .. ..  .. ..- attr(*, "predvars")= language list()
#>  .. .. .. .. ..  .. ..- attr(*, "dataClasses")= Named chr(0) 
#>  .. .. .. .. ..  .. .. ..- attr(*, "names")= chr(0) 
#>  .. .. .. .. .. $ xlev     : NULL
#>  .. .. .. .. .. $ contrasts: NULL
#>  .. ..  .. ..@ penalty   : NULL
#>  .. ..  .. ..@ formula   :Class 'formula'  language ~1
#>  .. .. .. .. ..  ..- attr(*, ".Environment")=<environment: 0x55e00116ca68> 
#>  .. @ response  : num [1:50] 1.354 1.574 2.105 3.17 0.946 ...
#>  .. @ n_obs     : int 50
#>  .. @ weights   : num [1:50] 1 1 1 1 1 1 1 1 1 1 ...
#>  .. @ offsets   :List of 2
#>  .. .. $ mu   : NULL
#>  .. .. $ sigma: NULL
#>  .. @ intercepts: Named logi [1:2] TRUE TRUE
#>  .. .. - attr(*, "names")= chr [1:2] "mu" "sigma"
#>  @ coefficients:List of 2
#>  .. $ mu   : num [1:2] 1.05 1.96
#>  .. $ sigma: num -1.43
#>  @ hyper       :List of 2
#>  .. $ mu   : list()
#>  .. $ sigma: list()
#>  @ loglik      : num 0.466
#>  @ objective   : num -0.466
#>  @ edf         :'data.frame':    2 obs. of  4 variables:
#>  .. $ parameter   : chr  "mu" "sigma"
#>  .. $ term        : chr  "linpar" "linpar"
#>  .. $ coefficients: int  2 1
#>  .. $ edf         : num  2 1
#>  @ fitted      :List of 2
#>  .. $ mu   : num [1:50] 1.57 1.78 2.17 2.83 1.44 ...
#>  .. $ sigma: num [1:50] 0.24 0.24 0.24 0.24 0.24 ...
#>  @ converged   : logi TRUE
#>  @ elapsed     : num 0.027
#>  @ history     :List of 3
#>  .. $ outer : NULL
#>  .. $ blocks:'data.frame':   1 obs. of  6 variables:
#>  ..  ..$ sweep     : int 1
#>  ..  ..$ block     : chr "smooth"
#>  ..  ..$ objective : num -0.466
#>  ..  ..$ change    : num 164
#>  ..  ..$ iterations: int 12
#>  ..  ..$ converged : logi TRUE
#>  .. $ inner :'data.frame':   11 obs. of  8 variables:
#>  ..  ..$ iteration: int [1:11] 1 2 3 4 5 6 7 8 9 10 ...
#>  ..  ..$ objective: num [1:11] 138.7 113.8 89.1 64.8 41.6 ...
#>  ..  ..$ score    : num [1:11] 3.709 0.999 0.996 0.99 0.972 ...
#>  ..  ..$ step     : num [1:11] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  ..$ rank     : int [1:11] 3 3 3 3 3 3 3 3 3 3 ...
#>  ..  ..$ route    : chr [1:11] "qr" "qr" "qr" "qr" ...
#>  ..  ..$ sweep    : int [1:11] 1 1 1 1 1 1 1 1 1 1 ...
#>  ..  ..$ block    : chr [1:11] "smooth" "smooth" "smooth" "smooth" ...
#>  @ methods     :List of 2
#>  .. $ smooth: <statmodels7::Iwls>
#>  ..  ..@ hessian      : chr "expected"
#>  ..  ..@ approx       : chr "bartlett"
#>  ..  ..@ decomposition: chr "qr"
#>  ..  ..@ maxit        : num 100
#>  ..  ..@ tol          : num 1e-06
#>  ..  ..@ step_halving : num 30
#>  .. $ sparse: chr(0) 
#>  @ call        : language statmod(formula = y ~ x, distrib = distributions7::gaussian1_distrib(),      data = sim2)
```
