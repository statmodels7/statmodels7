# What Kind of Block a Term Reports As

Returns the reading that a term receives in a summary: a structural
term's own parameters, the quantities of a break-point term, the
coefficients of an unpenalized term, a smooth's linear part and
smoothing parameter, a random effect's variance parameters, a
selection's survivors, or the coefficients of another penalized term.

## Usage

``` r
term_block_kind(term)
```

## Arguments

- term:

  A built term.

## Value

One of `"structural"`, `"breakpoint"`, `"parametric"`, `"smooth"`,
`"random"`, `"selection"`, `"penalized"`.

## Details

The classification is by the term's class and by its penalties, not by
its label, so a term given a name of its own is read the same way. A
structural term (a filter or a regime chain) and a break-point term are
classified first and report as `"structural"` and `"breakpoint"`
whatever penalties they carry. A random effect that carries a label is
read as `"random"`, because the covariance class that collects it
carries its penalty. For any other term the penalties are the ones that
the term declares through
[`modelterms7::term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.html):
a term with none, and with no labeled random effect in its subformulas,
is `"parametric"`; a term penalized over part of its parameters is read
as penalized, never as parametric, and as a selection when any of its
penalties has a kink.
