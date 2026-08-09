# statmodels7

[statmodels7](https://statmodels7.github.io/statmodels7/) installs and
attaches the packages of the
[statmodels7](https://statmodels7.github.io) toolkit, an S7 framework
for statistical modeling. Installing this package installs all of them,
and attaching it attaches all of them.

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/statmodels7")
```

Attaching it attaches the members, and says which arrived at which
version.

``` r

library(statmodels7)
#> -- Attaching the statmodels7 toolkit 0.3.0
#> v basis7         0.4.0    v numericals7    0.7.0
#> v distributions7 0.1.0    v optimizers7    0.1.0
#> v linkfunctions7 0.1.0    v parameters7    0.8.0
#> v modelterms7    0.13.0   v penalties7     0.5.0
```

To install or update every member afterwards:

``` r

statmodels7_update("install")
```

## What the toolkit is

Almost every R modeling package carries its own distributions and link
functions, written as internal helpers: a `switch` on a character
string, closures private to the package that owns them. They are not
objects, so nothing outside can reuse them, extend them, or ask them for
anything their author did not happen to need — and since each package
writes only what it needs, that is usually the density and the score,
sometimes the Hessian, and nothing beyond.

A Gamma is a fixed mathematical object. It should be written once,
correctly, with everything a modeling routine could want already
computed, and then everyone builds on top. The same argument applies to
an optimizer, whose stopping rule decides what a run means by *finished*
and is usually a number buried three levels down inside the function
that needs it.

``` r

statmodels7_packages()
#> [1] "basis7"         "distributions7" "linkfunctions7" "modelterms7"   
#> [5] "numericals7"    "optimizers7"    "parameters7"    "penalties7"
```

| package | what it provides |
|----|----|
| [linkfunctions7](https://statmodels7.github.io/linkfunctions7/) | link functions as objects, with exact derivatives to fourth order in both directions |
| [distributions7](https://statmodels7.github.io/distributions7/) | univariate and multivariate distributions with exact score, information and higher derivatives |
| [optimizers7](https://statmodels7.github.io/optimizers7/) | optimization algorithms as objects, with composable stopping rules |
| [basis7](https://statmodels7.github.io/basis7/) | basis expansions with derivatives, anchored integrals and exact Gram matrices |
| [parameters7](https://statmodels7.github.io/parameters7/) | constrained parameters as maps from an unconstrained vector, exact to fourth order |

## What this package does

Four functions, and an attach hook.

``` r

statmodels7_versions()
#>          package version
#> 1         basis7   0.4.0
#> 2 distributions7   0.1.0
#> 3 linkfunctions7   0.1.0
#> 4    modelterms7  0.13.0
#> 5    numericals7   0.7.0
#> 6    optimizers7   0.1.0
#> 7    parameters7   0.8.0
#> 8     penalties7   0.5.0
```

``` r

statmodels7_conflicts()
#> named list()
```

The members sit in `Imports` rather than `Depends`, which is what lets
the attaching be reported rather than merely happen: `Depends` would
attach them in the order the field lists them, with no message and no
way for a caller to see what arrived at which version.

## Related

The mathematical companion to the whole toolkit is the
[book](https://statmodels7.github.io/book/), which gives every formula
the packages implement with its derivation and its citation.
