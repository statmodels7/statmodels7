# Reject an Offset Buried Inside a Term

Stops where [`offset()`](https://rdrr.io/r/stats/offset.html) appears
anywhere other than as a term of the equation itself.

## Usage

``` r
reject_nested_offsets(eq, param)
```

## Arguments

- eq:

  A one-sided formula, with the equation's own offsets already taken out
  by
  [`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md).

- param:

  The parameter to which the equation belongs, named in the message.

## Value

`NULL`, invisibly, when the walk finds no offset. Signals an error
naming `param` and the offending call otherwise.

## Details

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
takes the top-level additive terms, which is where R recognizes an
offset. One written inside another term's formula, as in
`ridge(~ z + offset(o))`, `random(~ 1 + offset(o) | g)` or
`nl(a ~ 0 + ridge(~ g + offset(o)))`, reaches that term's own
`model.matrix` and is dropped there, as it would be dropped at the
equation level if it were left in. The block then carries the columns of
the model without the offset.

Such a term is rejected and is not moved to the equation, because the
meaning of an offset depends on where it is written. In a penalized
term's formula an offset would be a contribution to the equation's
predictor, and writing it at the equation level already expresses that.
In a subformula it would be a contribution to that parameter's own
chart, a different quantity. Picking one reading would fit the wrong
model in some cases, so the error names the place where the offset
belongs.

The whole expression tree is walked, so a term written later is covered.

## See also

[`split_offsets()`](https://statmodels7.github.io/statmodels7/reference/split_offsets.md)
