# Matching of numeric values

The `matchValues` method matches elements from `query` with those in
`target` using different matching approaches depending on parameter
`param`. Generally, `query` is expected to contain MS experimental
values (m/z and possibly retention time) while `target` reference
values. `query` and `target` can be `numeric`, a two dimensional array
(such as a `data.frame`, `matrix` or `DataFrame`), a
`SummarizedExperiment` or a `QFeatures`, `target` can in addition be a
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object. For `SummarizedExperiment`, the information for the matching is
expected to be in the object's `rowData`. For `QFeatures` matching is
performed for values present in the `rowData` of one of the object's
assays (which needs to be specified with the `assayQuery` parameter - if
a `QFeatures` is used as `target` the name of the assay needs to be
specified with parameter `assayTarget`). If `target` is a `Spectra`
matching is performed against spectra variables of this object and the
respective variable names need to be specified e.g. with `mzColname`
and/or `rtColname`. `matchMz` is an alias for `matchValues` to allow
backward compatibility.

Available `param` objects and corresponding matching approaches are:

- `ValueParam`: generic matching between values in `query` and `target`
  given acceptable differences expressed in `ppm` and `tolerance`. If
  `query` or `target` are not numeric, parameter `valueColname` has to
  be used to specify the name of the column that contains the values to
  be matched. The function returns a
  [`Matched()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
  object.

- `MzParam`: match query m/z values against reference compounds for
  which also m/z are known. Matching is performed similarly to the
  `ValueParam` above. If `query` or `target` are not numeric, the column
  name containing the values to be compared must be defined with
  `matchValues`' parameter `mzColname`, which defaults to `"mz"`.
  `MzParam` parameters `tolerance` and `ppm` allow to define the maximal
  acceptable (constant or m/z relative) difference between query and
  target m/z values.

- `MzRtParam`: match m/z **and** retention time values between `query`
  and `target`. Parameters `mzColname` and `rtColname` of the
  `matchValues` function allow to define the columns in `query` and
  `target` containing these values (defaulting to `c("mz", "mz")` and
  `c("rt", "rt")`, respectively). `MzRtParam` parameters `tolerance` and
  `ppm` have the same meaning as in `MzParam`; `MzRtParam` parameter
  `toleranceRt` allows to specify the maximal acceptable difference
  between query and target retention time values.

- `Mass2MzParam`: match m/z values against reference compounds for which
  only the (exact) mass is known. Before matching, m/z values are
  calculated from the compounds masses in the *target* table using the
  adducts specified via `Mass2MzParam` `adducts` parameter (defaults to
  `adducts = "[M+H]+"`). After conversion of adduct masses to m/z
  values, matching is performed similarly to `MzParam` (i.e. the same
  parameters `ppm` and `tolerance` can be used). If `query` is not
  `numeric`, parameter `mzColname` of `matchValues` can be used to
  specify the column containing the query's m/z values (defaults to
  `"mz"`). If `target` is a is not `numeric`, parameter `massColname`
  can be used to define the column containing the reference compound's
  masses (defaults to `"exactmass"`).

- `Mass2MzRtParam`: match m/z **and** retention time values against
  reference compounds for which the (exact) mass **and** retention time
  are known. Before matching, exact masses in `target` are converted to
  m/z values as for `Mass2MzParam`. Matching is then performed similarly
  to `MzRtParam`, i.e. m/z and retention times of entities are compared.
  With `matchValues`' parameters `mzColname`, `rtColname` and
  `massColname` the columns containing m/z values (in `query`),
  retention time values (in `query` and `target`) and exact masses (in
  `target`) can be specified.

- `Mz2MassParam`: input values for `query` and `target` are expected to
  be m/z values but matching is performed on exact masses calculated
  from these (based on the provided adduct definitions). In detail, m/z
  values in `query` are first converted to masses with the
  [`MetaboCoreUtils::mz2mass()`](https://rdrr.io/pkg/MetaboCoreUtils/man/mz2mass.html)
  function based on the adducts defined with `queryAdducts` (defaults to
  `"[M+H]+"`). The same is done for m/z values in `target` (adducts can
  be defined with `targetAdducts` which defaults to
  `"[M-H-]"). Matching is then performed on these converted values similarly to `ValueParam`. If `query`or`target`are not numeric, the column containing the m/z values can be specified with`matchValues`' parameter `mzColname`(defaults to`"mz"\`).

- `Mz2MassRtParam`: same as `Mz2MassParam` but with additional
  comparison of retention times between `query` and `target`. Parameters
  `rtColname` and `mzColname` of `matchValues` allow to specify which
  columns contain the retention times and m/z values, respectively.

## Usage

``` r
ValueParam(tolerance = 0, ppm = 5)

MzParam(tolerance = 0, ppm = 5)

Mass2MzParam(adducts = c("[M+H]+"), tolerance = 0, ppm = 5)

Mass2MzRtParam(adducts = c("[M+H]+"), tolerance = 0, ppm = 5, toleranceRt = 0)

MzRtParam(tolerance = 0, ppm = 0, toleranceRt = 0)

Mz2MassParam(
  queryAdducts = c("[M+H]+"),
  targetAdducts = c("[M-H]-"),
  tolerance = 0,
  ppm = 5
)

Mz2MassRtParam(
  queryAdducts = c("[M+H]+"),
  targetAdducts = c("[M+H]+"),
  tolerance = 0,
  ppm = 5,
  toleranceRt = 0
)

matchValues(query, target, param, ...)

# S4 method for class 'numeric,numeric,ValueParam'
matchValues(query, target, param)

# S4 method for class 'numeric,data.frameOrSimilar,ValueParam'
matchValues(
  query,
  target,
  param,
  valueColname = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,numeric,ValueParam'
matchValues(
  query,
  target,
  param,
  valueColname = character(),
  queryAssay = character()
)

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,ValueParam'
matchValues(
  query,
  target,
  param,
  valueColname = character(),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'numeric,numeric,Mass2MzParam'
matchValues(query, target, param)

# S4 method for class 'numeric,data.frameOrSimilar,Mass2MzParam'
matchValues(
  query,
  target,
  param,
  massColname = "exactmass",
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,numeric,Mass2MzParam'
matchValues(query, target, param, mzColname = "mz", queryAssay = character())

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,Mass2MzParam'
matchValues(
  query,
  target,
  param,
  mzColname = "mz",
  massColname = "exactmass",
  queryAssay = character(0),
  targetAssay = character(0)
)

# S4 method for class 'numeric,data.frameOrSimilar,MzParam'
matchValues(query, target, param, mzColname = "mz", targetAssay = character())

# S4 method for class 'numeric,Spectra,MzParam'
matchValues(query, target, param, mzColname = "mz", targetAssay = character())

# S4 method for class 'data.frameOrSimilar,numeric,MzParam'
matchValues(query, target, param, mzColname = "mz", queryAssay = character())

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,MzParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,Spectra,MzParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,Mass2MzRtParam'
matchValues(
  query,
  target,
  param,
  massColname = "exactmass",
  mzColname = "mz",
  rtColname = c("rt", "rt"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,MzRtParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  rtColname = c("rt", "rt"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,Spectra,MzRtParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  rtColname = c("rt", "rt"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'numeric,numeric,Mz2MassParam'
matchValues(query, target, param)

# S4 method for class 'numeric,data.frameOrSimilar,Mz2MassParam'
matchValues(query, target, param, mzColname = "mz", targetAssay = character())

# S4 method for class 'data.frameOrSimilar,numeric,Mz2MassParam'
matchValues(query, target, param, mzColname = "mz", queryAssay = character())

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,Mz2MassParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  queryAssay = character(),
  targetAssay = character()
)

# S4 method for class 'data.frameOrSimilar,data.frameOrSimilar,Mz2MassRtParam'
matchValues(
  query,
  target,
  param,
  mzColname = c("mz", "mz"),
  rtColname = c("rt", "rt"),
  queryAssay = character(),
  targetAssay = character()
)
```

## Arguments

- tolerance:

  for any `param` object: `numeric(1)` defining the maximal acceptable
  absolute difference in m/z (or in mass for `Mz2MassParam`) to consider
  them *matching*.

- ppm:

  for any `param` object: `numeric(1)` defining the maximal acceptable
  m/z-dependent (or mass-dependent for `Mz2MassParam`) difference (in
  parts-per-million) in m/z values to consider them to be *matching*.

- adducts:

  for `Mass2MzParam` or `Mass2MzRtParam`: either `character` with adduct
  names from
  [`MetaboCoreUtils::adducts()`](https://rdrr.io/pkg/MetaboCoreUtils/man/adductNames.html)
  or `data.frame` with a custom adduct definition. This parameter is
  used to calculate m/z from target compounds' masses. Custom adduct
  definitions can be passed to the adduct parameter in form of a
  `data.frame`. This `data.frame` is expected to have columns
  `"mass_add"` and `"mass_multi"` defining the *additive* and
  *multiplicative* part of the calculation. See
  [`MetaboCoreUtils::adducts()`](https://rdrr.io/pkg/MetaboCoreUtils/man/adductNames.html)
  for the expected format or use
  `MetaboCoreUtils::adductNames("positive")` and
  `MetaboCoreUtils::adductNames("negative")` for valid adduct names.

- toleranceRt:

  for `Mass2MzRtParam` or `MzRtParam`: `numeric(1)` defining the maximal
  acceptable absolute difference in retention time values to consider
  them them *matching*.

- queryAdducts:

  for `Mz2MassParam`. Adducts used to derive mass values from query m/z
  values. The expected format is the same as that for parameter
  `adducts`.

- targetAdducts:

  for `Mz2MassParam`. Adducts used to derive mass values from target m/z
  values. The expected format is the same as that for parameter
  `adducts`.

- query:

  feature table containing information on MS1 features. Can be a
  `numeric`, `data.frame`, `DataFrame`, `matrix`, `SummarizedExperiment`
  or `QFeatures`. It is expected to contain m/z values and can contain
  also other variables. Matchings based on both m/z and retention time
  can be performed when a column with retention times is present in both
  `query` and `target`.

- target:

  compound table with metabolites to compare against. The expected types
  are the same as those for `query`.

- param:

  parameter object defining the matching approach and containing the
  settings for that approach. See description above for details.

- ...:

  currently ignored.

- valueColname:

  `character` specifying the name of the column in `query` or/and the
  one in `target`with the desired values for the matching. This
  parameter should only be used when `param` is `valueParam` and in this
  case it must be provided (unless both `query` and `target` are
  `numeric`). It can be `character(1)` or `character(2)` in a similar
  way to `mzColname`.

- targetAssay:

  `character(1)` specifying the name of the assay of the provided
  `QFeatures` that should be used for the matching (values from this
  assay's `rowData` will be used for matching). Only used if `target` is
  an instance of a `QFeatures` object.

- queryAssay:

  `character(1)` specifying the name of the assay of the provided
  `QFeatures` that should be used for the matching (values from this
  assay's `rowData` will be used for matching). Only used if `query` is
  an instance of a `QFeatures` object.

- massColname:

  `character(1)` with the name of the column in `target` containing the
  mass of compounds. To be used when `param` is `Mass2MzParam` or
  `Mass2MzRtParam` (and target is not already `numeric` with the
  masses). Defaults to `massColname = "exactmass"`.

- mzColname:

  `character` specifying the name(s) of the column(s) in `query` or/and
  `target`with the m/z values. If one among `query` and `target` is
  `numeric` (and therefore there is no need to specify the column name)
  or `query` is not `numeric` and `param` is `Mass2MzParam` or
  `Mass2MzRtParam` (and therefore the name of the column with m/z needs
  only to be specified for `query`) then `mzColname` is expected to be
  `character(1)`. If both `query` and `target` are not numeric
  `mzColname` is expected to be `character(2)` (or `character(1)` and in
  this last case the two column names are assumed to be the same). If
  not specified the assumed default name for columns with m/z values is
  `"mz"`. If `target` is a
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object, the name of the spectra variable that should be used for the
  matching needs to be specified with `mzColname`.

- rtColname:

  `character(2)` with the name of the column containing the compounds
  retention times in `query` and the name for the one in `target`. It
  can also be `character(1)` if the two names are the same. To be used
  when `param` is `MzRtParam` or `Mass2MzRtParam`. Defaults to
  `rtColname = c("rt", "rt")`. If `target` is a
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object, the name of the spectra variable that should be used for the
  matching needs to be specified with `mzColname`.

## Value

[Matched](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
object representing the result.

Depending on the `param` object different *scores* representing the
quality of the match are provided. This comprises absolute as well as
relative differences (column/variables `"score"` and `"ppm_error"`
respectively). If `param` is a `Mz2MassParam`, `"score"` and
`"ppm_error"` represent differences of the compared masses (calculated
from the provided m/z values). If `param` an `MzParam`, `MzRtParam`,
`Mass2MzParam` or `Mass2MzRtParam`, `"score"` and `"ppm_error"`
represent absolute and relative differences of m/z values. Additionally,
if `param` is either an `MzRtParam` or `Mass2MzRtParam` differences
between query and target retention times for each matched element is
available in the column/variable `"score_rt"` in the returned `Matched`
object. Negative values of `"score"` (or `"score_rt"`) indicate that the
m/z or mass (or retention time) of the query element is smaller than
that of the target element.

## See also

[matchSpectra](https://rformassspectrometry.github.io/MetaboAnnotation/reference/matchSpectra.md)
or
[`CompareSpectraParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompareSpectraParam.md)
for spectra data matching

## Author

Andrea Vicini, Michael Witting

## Examples

``` r

library(MetaboCoreUtils)
#> 
#> Attaching package: ‘MetaboCoreUtils’
#> The following object is masked from ‘package:CompoundDb’:
#> 
#>     mass2mz
## Create a simple "target/reference" compound table
target_df <- data.frame(
   name = c("Tryptophan", "Leucine", "Isoleucine"),
   formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
   exactmass = c(204.089878, 131.094629, 131.094629)
)

## Create a "feature" table with m/z of features. We calculate m/z for
## certain adducts of some of the compounds in the reference table.
fts <- data.frame(
    feature_id = c("FT001", "FT002", "FT003"),
    mz = c(mass2mz(204.089878, "[M+H]+"),
           mass2mz(131.094629, "[M+H]+"),
           mass2mz(204.089878, "[M+Na]+") + 1e-6))

## Define the parameters for the matching
parm <- Mass2MzParam(
    adducts = c("[M+H]+", "[M+Na]+"),
    tolerance = 0,
    ppm = 20)
res <- matchValues(fts, target_df, parm)
res
#> Object of class Matched 
#> Total number of matches: 4 
#> Number of query objects: 3 (3 matched)
#> Number of target objects: 3 (3 matched)

## List the available variables/columns
colnames(res)
#> [1] "feature_id"       "mz"               "target_name"      "target_formula"  
#> [5] "target_exactmass" "adduct"           "score"            "ppm_error"       

## feature_id and mz are from the query data frame, while target_name,
## target_formula and target_exactmass are from the query object (columns
## from the target object have a prefix *target_* added to the original
## column name. Columns adduct, score and ppm_error represent the results
## of the matching: adduct the adduct/ion of the original compound for which
## the m/z matches, score the absolute difference of the query and target
## m/z and ppm_error the relative difference in m/z values.

## Get the full matching result:
matchedData(res)
#> DataFrame with 4 rows and 8 columns
#>      feature_id        mz target_name target_formula target_exactmass
#>     <character> <numeric> <character>    <character>        <numeric>
#> 1         FT001   205.097  Tryptophan     C11H12N2O2          204.090
#> 2         FT002   132.102     Leucine       C6H13NO2          131.095
#> 2.1       FT002   132.102  Isoleucine       C6H13NO2          131.095
#> 3         FT003   227.079  Tryptophan     C11H12N2O2          204.090
#>          adduct     score  ppm_error
#>     <character> <numeric>  <numeric>
#> 1        [M+H]+     0e+00 0.00000000
#> 2        [M+H]+     0e+00 0.00000000
#> 2.1      [M+H]+     0e+00 0.00000000
#> 3       [M+Na]+     1e-06 0.00440375

## We have thus matches of FT002 to two different compounds (but with the
## same mass).

## Individual columns can also be accessed with the $ operator:
res$feature_id
#> [1] "FT001" "FT002" "FT002" "FT003"
res$target_name
#> [1] "Tryptophan" "Leucine"    "Isoleucine" "Tryptophan"
res$ppm_error
#> [1] 0.000000000 0.000000000 0.000000000 0.004403752


## We repeat the matching requiring an exact match
parm <- Mass2MzParam(
    adducts = c("[M+H]+", "[M+Na]+"),
    tolerance = 0,
    ppm = 0)
res <- matchValues(fts, target_df, parm)
res
#> Object of class Matched 
#> Total number of matches: 3 
#> Number of query objects: 3 (2 matched)
#> Number of target objects: 3 (3 matched)

matchedData(res)
#> DataFrame with 4 rows and 8 columns
#>      feature_id        mz target_name target_formula target_exactmass
#>     <character> <numeric> <character>    <character>        <numeric>
#> 1         FT001   205.097  Tryptophan     C11H12N2O2          204.090
#> 2         FT002   132.102     Leucine       C6H13NO2          131.095
#> 2.1       FT002   132.102  Isoleucine       C6H13NO2          131.095
#> 3         FT003   227.079          NA             NA               NA
#>          adduct     score ppm_error
#>     <character> <numeric> <numeric>
#> 1        [M+H]+         0         0
#> 2        [M+H]+         0         0
#> 2.1      [M+H]+         0         0
#> 3            NA        NA        NA

## The last feature could thus not be matched to any compound.

## At last we use also different adduct definitions.
parm <- Mass2MzParam(
    adducts = c("[M+K]+", "[M+Li]+"),
    tolerance = 0,
    ppm = 20)
res <- matchValues(fts, target_df, parm)
res
#> Object of class Matched 
#> Total number of matches: 0 
#> Number of query objects: 3 (0 matched)
#> Number of target objects: 3 (0 matched)

matchedData(res)
#> DataFrame with 3 rows and 8 columns
#>    feature_id        mz target_name target_formula target_exactmass      adduct
#>   <character> <numeric> <character>    <character>        <numeric> <character>
#> 1       FT001   205.097          NA             NA               NA          NA
#> 2       FT002   132.102          NA             NA               NA          NA
#> 3       FT003   227.079          NA             NA               NA          NA
#>       score ppm_error
#>   <numeric> <numeric>
#> 1        NA        NA
#> 2        NA        NA
#> 3        NA        NA

## No matches were found.

## We can also match a "feature" table with a target data.frame taking into
## account both m/z and retention time values.
target_df <- data.frame(
  name = c("Tryptophan", "Leucine", "Isoleucine"),
  formula = c("C11H12N2O2", "C6H13NO2", "C6H13NO2"),
  exactmass = c(204.089878, 131.094629, 131.094629),
  rt = c(150, 140, 140)
)

fts <- data.frame(
  feature_id = c("FT001", "FT002", "FT003"),
  mz = c(mass2mz(204.089878, "[M+H]+"),
         mass2mz(131.094629, "[M+H]+"),
         mass2mz(204.089878, "[M+Na]+") + 1e-6),
  rt = c(150, 140, 150.1)
)

## Define the parameters for the matching
parm <- Mass2MzRtParam(
  adducts = c("[M+H]+", "[M+Na]+"),
  tolerance = 0,
  ppm = 20,
  toleranceRt = 0)

res <- matchValues(fts, target_df, parm)
res
#> Object of class Matched 
#> Total number of matches: 3 
#> Number of query objects: 3 (2 matched)
#> Number of target objects: 3 (3 matched)

## Get the full matching result:
matchedData(res)
#> DataFrame with 4 rows and 11 columns
#>      feature_id        mz        rt target_name target_formula target_exactmass
#>     <character> <numeric> <numeric> <character>    <character>        <numeric>
#> 1         FT001   205.097     150.0  Tryptophan     C11H12N2O2          204.090
#> 2         FT002   132.102     140.0     Leucine       C6H13NO2          131.095
#> 2.1       FT002   132.102     140.0  Isoleucine       C6H13NO2          131.095
#> 3         FT003   227.079     150.1          NA             NA               NA
#>     target_rt      adduct     score ppm_error  score_rt
#>     <numeric> <character> <numeric> <numeric> <numeric>
#> 1         150      [M+H]+         0         0         0
#> 2         140      [M+H]+         0         0         0
#> 2.1       140      [M+H]+         0         0         0
#> 3          NA          NA        NA        NA        NA

## FT003 could not be matched to any compound, FT002 was matched to two
## different compounds (but with the same mass).

## We repeat the matching allowing a positive tolerance for the matches
## between rt values

## Define the parameters for the matching
parm <- Mass2MzRtParam(
  adducts = c("[M+H]+", "[M+Na]+"),
  tolerance = 0,
  ppm = 20,
  toleranceRt = 0.1)

res <- matchValues(fts, target_df, parm)
res
#> Object of class Matched 
#> Total number of matches: 4 
#> Number of query objects: 3 (3 matched)
#> Number of target objects: 3 (3 matched)

## Get the full matching result:
matchedData(res)
#> DataFrame with 4 rows and 11 columns
#>      feature_id        mz        rt target_name target_formula target_exactmass
#>     <character> <numeric> <numeric> <character>    <character>        <numeric>
#> 1         FT001   205.097     150.0  Tryptophan     C11H12N2O2          204.090
#> 2         FT002   132.102     140.0     Leucine       C6H13NO2          131.095
#> 2.1       FT002   132.102     140.0  Isoleucine       C6H13NO2          131.095
#> 3         FT003   227.079     150.1  Tryptophan     C11H12N2O2          204.090
#>     target_rt      adduct     score  ppm_error  score_rt
#>     <numeric> <character> <numeric>  <numeric> <numeric>
#> 1         150      [M+H]+     0e+00 0.00000000       0.0
#> 2         140      [M+H]+     0e+00 0.00000000       0.0
#> 2.1       140      [M+H]+     0e+00 0.00000000       0.0
#> 3         150     [M+Na]+     1e-06 0.00440375       0.1

## Also FT003 was matched in this case

## It is also possible to match directly m/z values
mz1 <- c(12, 343, 23, 231)
mz2 <- mz1 + rnorm(4, sd = 0.001)

res <- matchValues(mz1, mz2, MzParam(tolerance = 0.001))

matchedData(res)
#> DataFrame with 4 rows and 4 columns
#>    query   target        score ppm_error
#>   <AsIs>   <AsIs>    <numeric> <numeric>
#> 1     12  11.9999  0.000108966  9.080580
#> 2    343 342.9999  0.000117242  0.341813
#> 3     23  23.0002 -0.000183083  7.960050
#> 4    231 231.0013 -0.001280555  5.543497

## Matching with a SummarizedExperiment or a QFeatures work analogously,
## only that the matching is performed on the object's `rowData`.

## Below we create a simple SummarizedExperiment with some random assay data.
## Note that results from a data preprocessing with the `xcms` package could
## be extracted as a `SummarizedExperiment` with the `quantify` method from
## the `xcms` package.
library(SummarizedExperiment)
se <- SummarizedExperiment(
    assays = matrix(rnorm(12), nrow = 3, ncol = 4,
    dimnames = list(NULL, c("A", "B", "C", "D"))),
    rowData = fts)

## We can now perform the matching of this SummarizedExperiment against the
## target_df as before.
res <- matchValues(se, target_df,
    param = Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
        tolerance = 0, ppm = 20))
res
#> Object of class Matched 
#> Total number of matches: 4 
#> Number of query objects: 3 (3 matched)
#> Number of target objects: 3 (3 matched)

## Getting the available columns
colnames(res)
#>  [1] "feature_id"       "mz"               "rt"               "target_name"     
#>  [5] "target_formula"   "target_exactmass" "target_rt"        "adduct"          
#>  [9] "score"            "ppm_error"       

## The query columns represent the columns of the object's `rowData`
rowData(se)
#> DataFrame with 3 rows and 3 columns
#>    feature_id        mz        rt
#>   <character> <numeric> <numeric>
#> 1       FT001   205.097     150.0
#> 2       FT002   132.102     140.0
#> 3       FT003   227.079     150.1

## matchedData also returns the query object's rowData along with the
## matching entries in the target object.
matchedData(res)
#> DataFrame with 4 rows and 10 columns
#>    feature_id        mz        rt target_name target_formula target_exactmass
#>   <character> <numeric> <numeric> <character>    <character>        <numeric>
#> 1       FT001   205.097     150.0  Tryptophan     C11H12N2O2          204.090
#> 2       FT002   132.102     140.0     Leucine       C6H13NO2          131.095
#> 3       FT002   132.102     140.0  Isoleucine       C6H13NO2          131.095
#> 4       FT003   227.079     150.1  Tryptophan     C11H12N2O2          204.090
#>   target_rt      adduct     score  ppm_error
#>   <numeric> <character> <numeric>  <numeric>
#> 1       150      [M+H]+     0e+00 0.00000000
#> 2       140      [M+H]+     0e+00 0.00000000
#> 3       140      [M+H]+     0e+00 0.00000000
#> 4       150     [M+Na]+     1e-06 0.00440375

## While `query` will return the full SummarizedExperiment.
query(res)
#> class: SummarizedExperiment 
#> dim: 3 4 
#> metadata(0):
#> assays(1): ''
#> rownames: NULL
#> rowData names(3): feature_id mz rt
#> colnames(4): A B C D
#> colData names(0):

## To illustrate use with a QFeatures object we first create a simple
## QFeatures object with two assays, `"ions"` representing the full feature
## data.frame and `"compounds"` a subset of it.
library(QFeatures)
#> Loading required package: MultiAssayExperiment
#> 
#> Attaching package: ‘QFeatures’
#> The following object is masked from ‘package:base’:
#> 
#>     sweep
qf <- QFeatures(list(ions = se, compounds = se[2,]))

## We can perform the same matching as before, but need to specify which of
## the assays in the QFeatures should be used for the matching. Below we
## perform the matching using the "ions" assay.
res <- matchValues(qf, target_df, queryAssay = "ions",
    param = Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
        tolerance = 0, ppm = 20))
res
#> Object of class Matched 
#> Total number of matches: 4 
#> Number of query objects: 3 (3 matched)
#> Number of target objects: 3 (3 matched)

## colnames returns now the colnames of the `rowData` of the `"ions"` assay.
colnames(res)
#>  [1] "feature_id"       "mz"               "rt"               "target_name"     
#>  [5] "target_formula"   "target_exactmass" "target_rt"        "adduct"          
#>  [9] "score"            "ppm_error"       

matchedData(res)
#> DataFrame with 4 rows and 10 columns
#>    feature_id        mz        rt target_name target_formula target_exactmass
#>   <character> <numeric> <numeric> <character>    <character>        <numeric>
#> 1       FT001   205.097     150.0  Tryptophan     C11H12N2O2          204.090
#> 2       FT002   132.102     140.0     Leucine       C6H13NO2          131.095
#> 3       FT002   132.102     140.0  Isoleucine       C6H13NO2          131.095
#> 4       FT003   227.079     150.1  Tryptophan     C11H12N2O2          204.090
#>   target_rt      adduct     score  ppm_error
#>   <numeric> <character> <numeric>  <numeric>
#> 1       150      [M+H]+     0e+00 0.00000000
#> 2       140      [M+H]+     0e+00 0.00000000
#> 3       140      [M+H]+     0e+00 0.00000000
#> 4       150     [M+Na]+     1e-06 0.00440375
```
