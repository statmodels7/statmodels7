# A Titled Rule for a Verbose Trace

One blank line, then a titled rule with the method that will do the work
named on the right, so that a reader of a running fit can see where one
step ends and the next begins.

## Usage

``` r
vb_rule(title, method = NULL, indent = 0L, char = "-")
```

## Arguments

- title:

  The step, e.g. `"outer 3"`.

- method:

  What runs it, or `NULL`.

- indent:

  How far in, in spaces.

- char:

  The rule's character.

## Value

Invisibly `NULL`; prints.

## Details

The trace has three nested things to say – which outer step, which pass
of the alternation inside it, and what each block did – and printed as
undifferentiated lines they are unreadable, as a panel fit with three
hyperparameters and 130 outer evaluations demonstrated. Naming the
method on every rule answers the question a reader of a slow fit
actually has, and that is the one running now.
