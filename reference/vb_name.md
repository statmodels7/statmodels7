# The Name of Whatever Is About to Run

Returns an optimizer's own name, a criterion's kind, or a short
description of a step that is not an optimizer.

## Usage

``` r
vb_name(x, default = "")
```

## Arguments

- x:

  An optimizer, a criterion, or `NULL`.

- default:

  The string returned when `x` is `NULL` or carries neither a name nor a
  kind.

## Value

A single string.
