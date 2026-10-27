# The Variables That a Term with a Moving Block Reads

Returns the data variables of a break-point term or of a nonlinear term:
the variable of the break-points and those of the subformulas for a
break-point term, and the data variables recorded when the term was
built for a nonlinear one.

## Usage

``` r
moving_term_vars(tm, name)
```

## Arguments

- tm:

  A
  [modelterms7::SegTerm](https://statmodels7.github.io/modelterms7/reference/SegTerm.html)
  or
  [modelterms7::NlTerm](https://statmodels7.github.io/modelterms7/reference/NlTerm.html).

- name:

  The name of the term in its equation.

## Value

A character vector.
