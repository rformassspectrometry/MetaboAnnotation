# Matching MS Spectra against a reference

`matchSpectra` compares experimental (*query*) MS2 spectra against
reference (*target*) MS2 spectra and reports matches with a similarity
that passing a specified threshold. The function performs the similarity
calculation between each query spectrum against each target spectrum.
Parameters `query` and `target` can be used to define the query and
target spectra, respectively, while parameter `param` allows to define
and configure the similarity calculation and matching condition.
Parameter `query` takes a
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html) object
while `target` can be either a
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html) object,
a [CompoundDb::CompDb](https://rdrr.io/pkg/CompoundDb/man/CompDb.html)
(reference library) object defined in the `CompoundDb` package or a
[CompAnnotationSource](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompAnnotationSource.md)
(e.g. a
[`CompDbSource()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompDbSource.md))
with the reference or connection information to a supported annotation
resource).

Some notes on performance and information on parallel processing are
provided in the vignette.

Currently supported parameter objects defining the matching are:

- `CompareSpectraParam`: the *generic* parameter object allowing to set
  all settings for the
  [`Spectra::compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  call that is used to perform the similarity calculation. This includes
  `MAPFUN` and `FUN` defining the peak-mapping and similarity
  calculation functions and `ppm` and `tolerance` to define an
  acceptable difference between m/z values of the compared peaks.
  Parameter `matchedPeaksCount` is also passed to `compareSpectra()`
  and, if set to `TRUE` (default is `FALSE`) will report the number of
  peaks defined to be *matching* by the `MAPFUN`. Additional parameters
  to the `compareSpectra` call can be passed along with `...`. See the
  help of
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for more information on these parameters. Importantly, if *msentropy*
  or a GNPS-like similarity calculation is used, `MAPFUN` should be
  selected accordingly (see section *Using alternative spectra
  similarity functions* in the package vignette for more information).
  By default, parameters `ppm` and `tolerance` are passed to the
  similarity calculation function, but if this function uses different
  parameters (e.g., `msentropy_similarity()` uses `ms2_tolerance_in_ppm`
  instead of `ppm`), these should be submitted to the
  `CompareSpectraParam()` function throught the `...` parameter.
  Parameters `requirePrecursor` (default `TRUE`) and
  `requirePrecursorPeak` (default `FALSE`) allow to pre-filter the
  target spectra prior to the actual similarity calculation for each
  individual query spectrum. Parameters `ppm` and `tolerance` are also
  used to define the maximal acceptable difference in precursor m/z if
  `requirePrecursor` or `requirePrecursorPeak` are set to `TRUE`. Target
  spectra can also be pre-filtered based on retention time if parameter
  `toleranceRt` is set to a value different than the default
  `toleranceRt = Inf`. Only target spectra with a retention time within
  the query's retention time +/- (`toleranceRt` + `percentRt`% of the
  query's retention time) are considered. Note that while for `ppm` and
  `tolerance` only a single value is accepted, `toleranceRt` and
  `percentRt` can be also of length equal to the number of query spectra
  hence allowing to define different rt boundaries for each query
  spectrum. While these pre-filters can considerably improve
  performance, it should be noted that no matches will be found between
  query and target spectra with missing values in the considered
  variable (precursor m/z or retention time). For target spectra without
  retention times (such as for `Spectra` from a public reference
  database such as MassBank) the default `toleranceRt = Inf` should thus
  be used. Finally, parameter `THRESHFUN` allows to define a function to
  be applied to the similarity scores to define which matches to report.
  See below for more details.

- `MatchForwardReverseParam`: performs spectra matching as with
  `CompareSpectraParam` but reports, similar to MS-DIAL, also the
  *reverse* similarity score and the *presence ratio*. Please refer to
  the documentation of `CompareSpectraParam` for explanation of the
  parameters. With `MatchForwardReverseParam`, the matching of query
  spectra to target spectra is performed by considering all peaks from
  the query and all peaks from the target (reference) spectrum (i.e.
  *forward* matching using an *outer join*-based peak matching
  strategy). For matching spectra also the *reverse* similarity is
  calculated considering only peaks present in the target (reference)
  spectrum (i.e. using a *right join*-based peak matching). This is
  reported as spectra variable `"reverse_score"`. In addition, the ratio
  between the number of matched peaks and the total number of peaks in
  the target (reference) spectra is reported as the *presence ratio*
  (spectra variable `"presence_ratio"`) and the total number of matched
  peaks as `"matched_peaks_count"`. See examples below for details.
  Parameter `THRESHFUN_REVERSE` allows to define an additional
  *threshold function* to filter matches. If `THRESHFUN_REVERSE` is
  defined only matches with a spectra similarity fulfilling both
  `THRESHFUN` **and** `THRESHFUN_REVERSE` are returned. With the default
  `THRESHFUN_REVERSE = NULL` all matches passing `THRESHFUN` are
  reported.

## Usage

``` r
# S4 method for class 'Spectra,CompDbSource,Param'
matchSpectra(
  query,
  target,
  param,
  BPPARAM = BiocParallel::SerialParam(),
  addOriginalQueryIndex = TRUE
)

CompareSpectraParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 5,
  FUN = MsCoreUtils::ndotproduct,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.7),
  toleranceRt = Inf,
  percentRt = 0,
  matchedPeaksCount = FALSE,
  ...
)

MatchForwardReverseParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 5,
  FUN = MsCoreUtils::ndotproduct,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.7),
  THRESHFUN_REVERSE = NULL,
  toleranceRt = Inf,
  percentRt = 0,
  ...
)

# S4 method for class 'Spectra,Spectra,CompareSpectraParam'
matchSpectra(
  query,
  target,
  param,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SerialParam(),
  addOriginalQueryIndex = TRUE
)

# S4 method for class 'Spectra,CompDb,Param'
matchSpectra(
  query,
  target,
  param,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SerialParam(),
  addOriginalQueryIndex = TRUE
)

# S4 method for class 'Spectra,Spectra,MatchForwardReverseParam'
matchSpectra(
  query,
  target,
  param,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SerialParam(),
  addOriginalQueryIndex = TRUE
)
```

## Arguments

- query:

  for `matchSpectra`:
  [Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object with the query spectra.

- target:

  for `matchSpectra`:
  [Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html),
  [CompoundDb::CompDb](https://rdrr.io/pkg/CompoundDb/man/CompDb.html)
  or object extending
  [CompAnnotationSource](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompAnnotationSource.md)
  (such as
  [CompDbSource](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompDbSource.md))
  with the target (reference) spectra to compare `query` against.

- param:

  for `matchSpectra`: parameter object (such as `CompareSpectraParam`)
  defining the settings for the matching.

- BPPARAM:

  for `matchSpectra`: parallel processing setup (see the `BiocParallel`
  package for more information). Parallel processing is disabled by
  default (with the default setting `BPPARAM = SerialParam()`).

- addOriginalQueryIndex:

  for
  [`matchSpectra()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/matchSpectra.md):
  `logical(1)` whether an additional spectra variable
  `".original_query_index"` should be added to the `query` `Spectra`
  object providing the index of the spectrum in this originally provided
  object. This spectra variable can be useful to link back to the
  original `Spectra` object if the `MatchedSpectra` object gets
  subsetted/processed.

- MAPFUN:

  `function` used to map peaks between the compared spectra. Defaults
  for `CompareSpectraParam` to
  [`Spectra::joinPeaks()`](https://rdrr.io/pkg/Spectra/man/joinPeaks.html).
  See
  [`Spectra::compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  for details.

- tolerance:

  `numeric(1)` for an absolute maximal accepted difference between m/z
  values. This will be used in `compareSpectra` as well as for eventual
  precursor m/z matching.

- ppm:

  `numeric(1)` for a relative, m/z-dependent, maximal accepted
  difference between m/z values. This will be used in `compareSpectra`
  as well as for eventual precursor m/z matching.

- FUN:

  `function` used to calculate similarity between spectra. Defaults for
  `CompareSpectraParam` to
  [`MsCoreUtils::ndotproduct()`](https://rdrr.io/pkg/MsCoreUtils/man/distance.html).
  See
  [`MsCoreUtils::ndotproduct()`](https://rdrr.io/pkg/MsCoreUtils/man/distance.html)
  for details.

- requirePrecursor:

  `logical(1)` whether only target spectra are considered in the
  similarity calculation with a precursor m/z that matches the precursor
  m/z of the query spectrum (considering also `ppm` and `tolerance`).
  With `requirePrecursor = TRUE` (the default) the function will
  complete much faster, but will not find any hits for target (or query
  spectra) with missing precursor m/z. It is suggested to check first
  the availability of the precursor m/z in `target` and `query`.

- requirePrecursorPeak:

  `logical(1)` whether only target spectra will be considered in the
  spectra similarity calculation that have a peak with an m/z matching
  the precursor m/z of the query spectrum. Defaults to
  `requirePrecursorPeak = FALSE`. It is suggested to check first the
  availability of the precursor m/z in `query`, as no match will be
  reported for query spectra with missing precursor m/z.

- THRESHFUN:

  `function` applied to the similarity score to define which target
  spectra are considered *matching*. Defaults to
  `THRESHFUN = function(x) which(x >= 0.7)` hence selects all target
  spectra matching a query spectrum with a similarity higher or equal
  than `0.7`. Any function that takes a numeric vector with similarity
  scores from the comparison of a query spectrum with all target spectra
  (as returned by
  [`Spectra::compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html))
  as input and returns a `logical` vector (same dimensions as the
  similarity scores) or an integer with the matches is supported.

- toleranceRt:

  `numeric` of length 1 or equal to the number of query spectra defining
  the maximal accepted (absolute) difference in retention time between
  query and target spectra. By default (with `toleranceRt = Inf`) the
  retention time-based filter is not considered. See help of
  `CompareSpectraParam` above for more information.

- percentRt:

  `numeric` of length 1 or equal to the number of query spectra defining
  the maximal accepted relative difference in retention time between
  query and target spectra expressed in percentage of the query rt. For
  `percentRt = 10`, similarities are defined between the query spectrum
  and all target spectra with a retention time within query rt +/- 10%
  of the query. By default (with `toleranceRt = Inf`) the retention
  time-based filter is not considered. Thus, to consider the `percentRt`
  parameter, `toleranceRt` should be set to a value different than that.
  See help of `CompareSpectraParam` above for more information.

- matchedPeaksCount:

  `logical(1)` for `CompareSpectraParam()`: whether also the number of
  matching peaks should be reported (in column `"matched_peaks_count"`).
  This number represents the number of peaks reported *matching* by the
  `MAPFUN`.

- ...:

  for `CompareSpectraParam`: additional parameters passed along to the
  [`Spectra::compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  call, including eventual additional parameters of the selected mapping
  or similarity calculation functions.

- THRESHFUN_REVERSE:

  for `MatchForwardReverseParam`: optional additional *thresholding
  function* to filter the results on the reverse score. If specified the
  same format than `THRESHFUN` is expected.

- rtColname:

  `character(2)` with the name of the spectra variable containing the
  retention time information for compounds to be used in retention time
  matching (only used if `toleranceRt` is not `Inf`). It can also be
  `character(1)` if the two names are the same. Defaults to
  `rtColname = c("rtime", "rtime")`.

## Value

`matchSpectra` returns a
[`MatchedSpectra()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/MatchedSpectra.md)
object with the matching results. If `target` is a
`CompAnnotationSource` only matching target spectra will be reported.

Constructor functions return an instance of the class.

## Author

Johannes Rainer, Michael Witting

## Examples

``` r

library(Spectra)
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:MetaboAnnotation’:
#> 
#>     endoapply
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: BiocParallel
library(MsDataHub)
## Load a test file from *MsDataHub*
fl <- MsDataHub::PestMix1_DDA.mzML()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)

## subset to selected spectra.
pest_ms2 <- pest_ms2[c(808, 809, 945:955)]

## Load a small example MassBank data set
load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))

## Match spectra with the default similarity score (normalized dot product)
csp <- CompareSpectraParam(requirePrecursor = TRUE, ppm = 10)
mtches <- matchSpectra(pest_ms2, minimb, csp)

mtches
#> Object of class MatchedSpectra 
#> Total number of matches: 16 
#> Number of query objects: 13 (5 matched)
#> Number of target objects: 100 (11 matched)

## Are there any matching spectra for the first query spectrum?
mtches[1]
#> Object of class MatchedSpectra 
#> Total number of matches: 0 
#> Number of query objects: 1 (0 matched)
#> Number of target objects: 100 (0 matched)
## No

## And for the second query spectrum?
mtches[2]
#> Object of class MatchedSpectra 
#> Total number of matches: 4 
#> Number of query objects: 1 (1 matched)
#> Number of target objects: 100 (4 matched)
## The second query spectrum matches 4 target spectra. The scores for these
## matches are:
mtches[2]$score
#> [1] 0.7869556 0.8855473 0.7234894 0.7219942

## To access the score for the full data set
mtches$score
#>  [1]        NA 0.7869556 0.8855473 0.7234894 0.7219942        NA 0.7769746
#>  [8] 0.7577286        NA 0.7433718 0.7019807 0.7081274        NA 0.7320465
#> [15] 0.8106258 0.7290458 0.8168876 0.7247800 0.7412586 0.7198787        NA
#> [22]        NA        NA        NA

## Below we use a THRESHFUN that returns for each query spectrum the (first)
## best matching target spectrum.
csp <- CompareSpectraParam(requirePrecursor = FALSE, ppm = 10,
    THRESHFUN = function(x) which.max(x))
mtches <- matchSpectra(pest_ms2, minimb, csp)
mtches
#> Object of class MatchedSpectra 
#> Total number of matches: 13 
#> Number of query objects: 13 (13 matched)
#> Number of target objects: 100 (10 matched)

## Each of the query spectra is matched to one target spectrum
length(mtches)
#> [1] 13
matches(mtches)
#>    query_idx target_idx        score
#> 1          1          1 0.000000e+00
#> 2          2         73 8.855473e-01
#> 3          3          2 6.313687e-01
#> 4          4         44 7.769746e-01
#> 5          5         74 1.772117e-05
#> 6          6          2 7.433718e-01
#> 7          7          5 1.906998e-03
#> 8          8         53 8.168876e-01
#> 9          9         44 7.412586e-01
#> 10        10         86 4.085289e-04
#> 11        11         53 4.323403e-01
#> 12        12         47 3.469648e-03
#> 13        13         71 7.612480e-06

## Match spectra considering also measured retention times. This requires
## that both query and target spectra have non-missing retention times.
rtime(pest_ms2)
#>  [1] 361.651 361.741 377.609 377.699 378.120 378.539 378.779 378.869 378.959
#> [10] 379.379 380.059 380.609 381.029
rtime(minimb)
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA

## Target spectra don't have retention times. Below we artificially set
## retention times to show how an additional retention time filter would
## work.
rtime(minimb) <- rep(361, length(minimb))

## Matching spectra requiring a matching precursor m/z and the difference
## of retention times between query and target spectra to be <= 2 seconds.
csp <- CompareSpectraParam(requirePrecursor = TRUE, ppm = 10,
    toleranceRt = 2)
mtches <- matchSpectra(pest_ms2, minimb, csp)
mtches
#> Object of class MatchedSpectra 
#> Total number of matches: 4 
#> Number of query objects: 13 (1 matched)
#> Number of target objects: 100 (4 matched)
matches(mtches)
#>   query_idx target_idx     score
#> 1         2         70 0.7869556
#> 2         2         73 0.8855473
#> 3         2         75 0.7234894
#> 4         2         76 0.7219942

## Note that parameter `rtColname` can be used to define different spectra
## variables with retention time information (such as retention indices etc).

## A `CompDb` compound annotation database could also be used with
## parameter `target`. Below we load the test `CompDb` database from the
## `CompoundDb` Bioconductor package.
library(CompoundDb)
#> Loading required package: AnnotationFilter
fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
cdb <- CompDb(fl)
res <- matchSpectra(pest_ms2, cdb, CompareSpectraParam())

## We do however not find any matches since the used compound annotation
## database contains only a very small subset of the MassBank.
res
#> Object of class MatchedSpectra 
#> Total number of matches: 0 
#> Number of query objects: 13 (0 matched)
#> Number of target objects: 70 (0 matched)

## As `target` we have now however the MS2 spectra data from the compound
## annotation database
target(res)
#> MSn data (Spectra) with 70 spectra in a MsBackendCompDb backend:
#>       msLevel precursorMz  polarity
#>     <integer>   <numeric> <integer>
#> 1           2      179.07         1
#> 2           2      179.07         1
#> 3           2      179.07         1
#> 4           2      179.07         1
#> 5           2      179.07         1
#> ...       ...         ...       ...
#> 66          2     337.091         1
#> 67          2     337.091         1
#> 68          2     337.091         1
#> 69          2     337.091         1
#> 70          2     337.091         1
#>  ... 47 more variables/columns.
#>  Use  'spectraVariables' to list all of them.
#>  data source: MassBank 
#>  version: 2020.09 
#>  organism: NA 

## See the package vignette for details, descriptions and more examples,
## also on how to retrieve e.g. MassBank reference databases from
## Bioconductor's AnnotationHub.
```
