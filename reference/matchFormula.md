# Chemical Formula Matching

The `matchFormula` method matches chemical formulas from different
inputs (parameter `query` and `target`). Before comparison all formulas
are normalized using
[`MetaboCoreUtils::standardizeFormula()`](https://rdrr.io/pkg/MetaboCoreUtils/man/standardizeFormula.html).
Inputs can be either a `character` or `data.frame` containing a column
with formulas. In case of `data.frame`s parameter `formulaColname` needs
to be used to specify the name of the column containing the chemical
formulas.

## Usage

``` r
matchFormula(query, target, ...)

# S4 method for class 'character,character'
matchFormula(query, target, BPPARAM = SerialParam())

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar'
matchFormula(
  query,
  target,
  formulaColname = c("formula", "formula"),
  BPPARAM = SerialParam()
)

# S4 method for class 'character,data.frameOrSimilar'
matchFormula(
  query,
  target,
  formulaColname = "formula",
  BPPARAM = SerialParam()
)

# S4 method for class 'data.frameOrSimilar,character'
matchFormula(
  query,
  target,
  formulaColname = "formula",
  BPPARAM = SerialParam()
)
```

## Arguments

- query:

  `character` or `data.frame` with chemical formulas to search.

- target:

  `character` or `data.frame` with chemical formulas to compare against.

- ...:

  currently ignored

- BPPARAM:

  parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

- formulaColname:

  `character` with the name of the column containing chemical formulas.
  Can be of length 1 if both `query` and `target` are `data.frame`s and
  the name of the column with chemical formulas is the same for both. If
  different columns are used, `formulaColname[1]` can be used to define
  the column name in `query` and `formulaColname[2]` the one of
  `target`.

## Value

[Matched](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
object representing the result.

## Author

Michael Witting

## Examples

``` r

## input formula
query <- c("H12C6O6", "C11H12O2", "HN3")
target <- c("HCl", "C2H4O", "C6H12O6")

query_df <- data.frame(
    formula = c("H12C6O6", "C11H12O2", "HN3"),
    name = c("A", "B", "C")
)
target_df <- data.frame(
    formula = c("HCl", "C2H4O", "C6H12O6"),
    name = c("D", "E", "F")
)

## character vs character
matches <- matchFormula(query, target)
matchedData(matches)
#> DataFrame with 3 rows and 3 columns
#>            query  target     score
#>           <AsIs>  <AsIs> <numeric>
#> C6H12O6  H12C6O6 C6H12O6         1
#> NA      C11H12O2      NA        NA
#> NA.1         HN3      NA        NA

## data.frame vs data.frame
matches <- matchFormula(query_df, target_df)
matchedData(matches)
#> DataFrame with 3 rows and 5 columns
#>       formula        name target_formula target_name     score
#>   <character> <character>    <character> <character> <numeric>
#> 1     H12C6O6           A        C6H12O6           F         1
#> 2    C11H12O2           B             NA          NA        NA
#> 3         HN3           C             NA          NA        NA
## data.frame vs character
matches <- matchFormula(query_df, target)
matchedData(matches)
#> DataFrame with 3 rows and 4 columns
#>       formula        name  target     score
#>   <character> <character>  <AsIs> <numeric>
#> 1     H12C6O6           A C6H12O6         1
#> 2    C11H12O2           B      NA        NA
#> 3         HN3           C      NA        NA
## character vs data.frame
matches <- matchFormula(query, target_df)
matchedData(matches)
#> DataFrame with 3 rows and 4 columns
#>         query target_formula target_name     score
#>        <AsIs>    <character> <character> <numeric>
#> 3     H12C6O6        C6H12O6           F         1
#> NA   C11H12O2             NA          NA        NA
#> NA.1      HN3             NA          NA        NA
```
