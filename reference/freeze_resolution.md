# The Outer Optimizer as the Fit Stores It

Returns the optimizer of the outer search with the resolution of its
line search stored as a number.

## Usage

``` r
freeze_resolution(optimizer, ours, value)
```

## Arguments

- optimizer:

  The optimizer that ran the search.

- ours:

  `TRUE` where
  [`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
  chose the optimizer and set the function, `FALSE` where the caller
  supplied it.

- value:

  The resolution that the function returns at the end of the search.

## Value

The optimizer.

## Details

[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md)
gives the line search of the optimizer it chooses a function that
returns the resolution of the criterion at the moment of each call, so
that the line search follows the readings taken during the run. The
function is created inside
[`outer_fit()`](https://statmodels7.github.io/statmodels7/reference/outer_fit.md),
and its environment is the frame of that call, which holds the design,
the inner fit and the caches of the derivatives. Stored in
`methods$search`, the function takes the frame into every serialization
of the fit. After the search no reading is taken, the function returns
the same number at every call, and the optimizer is stored with that
number. A value that is not finite is stored as zero, the resolution of
a line search with none.

An optimizer that the caller supplied is returned unchanged, and so is
one whose line search has no function in `resolution`.
