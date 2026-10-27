# Reject an Argument That a Method Does Not Take

Signals an error when a method that has `...` for compatibility with its
generic receives an argument in it. The message names each argument, as
in `vcov() takes no argument 'x'.`, and shows `<unnamed>` for an
argument without a name. The check keeps a misspelled or nonexistent
argument from being accepted and ignored.

## Usage

``` r
check_dots(fn, ...)
```

## Arguments

- fn:

  The name of the method's generic, for the message.

- ...:

  The dots of the method.

## Value

`NULL`, invisibly, when the dots are empty.
