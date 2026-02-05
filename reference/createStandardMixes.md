# Create Standard Mixes from a Matrix of Standard Compounds

The `createStandardMixes` function defines groups (mixes) of compounds
(standards) with dissimilar m/z values. The expected size of the groups
can be defined with parameters `max_nstd` and `min_nstd` and the minimum
required difference between m/z values within each group with parameter
`min_diff`. The group assignment will be reported in an additional
column in the result data frame.

## Usage

``` r
createStandardMixes(
  x,
  max_nstd = 10,
  min_nstd = 5,
  min_diff = 2,
  iterativeRandomization = FALSE
)
```

## Arguments

- x:

  `numeric` matrix with row names representing the compounds and columns
  representing different adducts. Such a matrix with m/z values for
  different adducts for compounds could e.g. be created with the
  [`MetaboCoreUtils::mass2mz()`](https://rdrr.io/pkg/MetaboCoreUtils/man/mass2mz.html)
  function.

- max_nstd:

  `numeric` number of maximum standards per group.

- min_nstd:

  `numeric` number of minimum standards per group. Only needed when
  using `iterativeRandomization = TRUE`.

- min_diff:

  `numeric` Minimum difference for considering two values as distinct.

- iterativeRandomization:

  `logical` default `FALSE`. If set to `TRUE`, `createStandardMixes`
  will randomly rearrange the rows of `x` until the user inputs are
  satisfied.

## Value

`data.frame` created by adding a column `group` to the input `x` matrix,
comprising the group number for each compound.

## Details

Users should be aware that because the function iterates through `x`,
the compounds at the bottom of the matrix are more complicated to group,
and there is a possibility that some compounds will not be grouped with
others. We advise specifyiong `iterativeRandomization = TRUE` even if it
takes more time.

## Author

Philippine Louail

## Examples

``` r

## Iterative grouping only
x <- matrix(c(135.0288, 157.0107, 184.0604, 206.0424, 265.1118, 287.0937,
              169.0356, 191.0176, 468.9809, 490.9628, 178.0532, 200.0352),
            ncol = 2, byrow = TRUE,
            dimnames = list(c("Malic Acid", "Pyridoxic Acid", "Thiamine",
                                "Uric acid", "dUTP", "N-Formyl-L-methionine"),
                             c("adduct_1", "adduct_2")))
result <- createStandardMixes(x, max_nstd = 3, min_diff = 2)

## Randomize grouping
set.seed(123)
x <- matrix(c(349.0544, 371.0363, 325.0431, 347.0251, 581.0416, 603.0235,
              167.0564, 189.0383, 150.0583, 172.0403, 171.0053, 192.9872,
              130.0863, 152.0682, 768.1225, 790.1044),
            ncol = 2, byrow = TRUE,
            dimnames = list(c("IMP", "UMP", "UDP-glucuronate",
                                "1-Methylxanthine", "Methionine",
                                "Dihydroxyacetone phosphate",
                                "Pipecolic acid", "CoA"),
                             c("[M+H]+", "[M+Na]+")))
result <- createStandardMixes(x, max_nstd = 4, min_nstd = 3, min_diff = 2,
                               iterativeRandomization = TRUE)
```
