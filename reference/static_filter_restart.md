# Set the Constants of a Statically Started Filter From Its Static Model

Reads the predictor of the static model at the first row of each group
of a statically started filter and stores it, with the model, on the
term; for a joint filter it does so for each member, on that member's
parameter.

## Usage

``` r
static_filter_restart(tm, param, fit_s, newdata)
```

## Arguments

- tm:

  A model term; any term other than a statically started filter is
  returned unchanged.

- param:

  The distribution parameter in whose equation the term sits.

- fit_s:

  The fitted static model.

- newdata:

  `NULL` for the fitting rows, or the rows on which the term was
  rebuilt.

## Value

The term.
