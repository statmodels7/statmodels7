# A Titled Rule for a Verbose Trace

One blank line, then a titled rule with the name of the method that does
the work on the right, so that a reader of a running fit can see where
one step ends and the next begins.

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

A verbose trace reports three nested levels: the outer step, the pass of
the alternation inside it, and the work of each block. The rule
separates these levels, and naming the method on every rule shows which
method runs at that step.
