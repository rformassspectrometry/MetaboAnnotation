# Representation of Spectra matches

Matches between query and target spectra can be represented by the
`MatchedSpectra` object. Functions like the
[`matchSpectra()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/matchSpectra.md)
function will return this type of object. By default, all data accessors
work as *left joins* between the *query* and the *target* spectra, i.e.
values are returned for each *query* spectrum with eventual duplicated
entries (values) if the query spectrum matches more than one target
spectrum.

## Usage

``` r
MatchedSpectra(
  query = Spectra(),
  target = Spectra(),
  matches = data.frame(query_idx = integer(), target_idx = integer(), score = numeric())
)

# S4 method for class 'MatchedSpectra'
spectraVariables(object)

# S4 method for class 'MatchedSpectra'
queryVariables(object)

# S4 method for class 'MatchedSpectra'
targetVariables(object)

# S4 method for class 'MatchedSpectra'
colnames(x)

# S4 method for class 'MatchedSpectra'
x$name

# S4 method for class 'MatchedSpectra'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MatchedSpectra'
matchedData(object, columns = spectraVariables(object), ...)

# S4 method for class 'MatchedSpectra'
addProcessing(object, FUN, ..., spectraVariables = character())

# S4 method for class 'MatchedSpectra'
plotSpectraMirror(
  x,
  xlab = "m/z",
  ylab = "intensity",
  main = "",
  scalePeaks = FALSE,
  ...
)

# S4 method for class 'MatchedSpectra,MsBackend'
setBackend(object, backend, ...)
```

## Arguments

- query:

  `Spectra` with the query spectra.

- target:

  `Spectra` with the spectra against which `query` has been matched.

- matches:

  `data.frame` with columns `"query_idx"` (`integer`), `"target_idx"`
  (`integer`) and `"score"` (`numeric`) representing the *n:m* mapping
  of elements between the `query` and the `target` `Spectra`.

- object:

  `MatchedSpectra` object.

- x:

  `MatchedSpectra` object.

- name:

  for `$`: the name of the spectra variable to extract.

- columns:

  for `spectraData`: `character` vector with spectra variable names that
  should be extracted.

- ...:

  for `addProcessing`: additional parameters for the function `FUN`. For
  `plotSpectraMirror`: additional parameters passed to the plotting
  functions.

- FUN:

  for `addProcessing`: function to be applied to the peak matrix of each
  spectrum in `object`. See
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for more details.

- spectraVariables:

  for `addProcessing`: `character` with additional spectra variables
  that should be passed along to the function defined with `FUN`. See
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for details.

- xlab:

  for `plotSpectraMirror`: the label for the x-axis.

- ylab:

  for `plotSpectraMirror`: the label for the y-axis.

- main:

  for `plotSpectraMirror`: an optional title for each plot.

- scalePeaks:

  for `plotSpectraMirror`: `logical(1)` if peak intensities (per
  spectrum) should be scaled to a total sum of one (per spectrum) prior
  to plotting.

- backend:

  for `setBackend`: instance of an object extending
  [Spectra::MsBackend](https://rdrr.io/pkg/Spectra/man/MsBackend.html).
  See help for
  [`Spectra::setBackend()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for more details.

## Value

See individual method desciption above for details.

## Creation, subset and filtering

`MatchedSpectra` objects are the result object from the
[`matchSpectra()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/matchSpectra.md).
While generally not needed, `MatchedSpectra` objects can also be created
with the `MatchedSpectra` function providing the `query` and `target`
`Spectra` objects as well as a `data.frame` with the *matches* between
query and target elements. This data frame is expected to have columns
`"query_idx"`, `"target_idx"` with the `integer` indices of query and
target objects that are *matched* and a column `"score"` with a
`numeric` score for the match.

`MatchedSpectra` objects can be subset using:

- `[` subset the `MatchedSpectra` selecting `query` spectra to keep with
  parameter `i`. The `target` spectra will by default be returned as-is.

- `pruneTarget` *cleans* the `MatchedSpectra` object by removing
  non-matched target spectra.

In addition, `MatchedSpectra` can be filtered with any of the filtering
approaches defined for
[`Matched()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
objects:
[`SelectMatchesParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md),
[`TopRankedMatchesParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
or
[`ScoreThresholdParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md).

## Extracting data

- `$` extracts a single spectra variable from the `MatchedSpectra` `x`.
  Use `spectraVariables` to get all available spectra variables. Prefix
  `"target_"` is used for spectra variables from the *target* `Spectra`.
  The matching scores are available as *spectra variable* `"score"`.
  Similar to a left join between the query and target spectra, this
  function returns a value for each query spectrum with eventual
  duplicated values for query spectra matching more than one target
  spectrum. If spectra variables from the target spectra are extracted,
  an `NA` is reported for *query* spectra that don't match any target
  spectra. See examples below for more details.

- `length` returns the number of **query** spectra.

- `matchedData` same as `spectraData` below.

- `query` returns the *query* `Spectra`.

- `queryVariables` returns the `spectraVariables` of *query*.

- `spectraData` returns spectra variables from the query and/or target
  `Spectra` as a `DataFrame`. Parameter `columns` allows to define which
  variables should be returned (defaults to
  `columns = spectraVariables(object)`), spectra variable names of the
  target spectra need to be prefixed with `target_` (e.g.
  `target_msLevel` to get the MS level from target spectra). The score
  from the matching function is returned as spectra variable `"score"`.
  Similar to `$`, this function performs a *left join* of spectra
  variables from the *query* and *target* spectra returning all values
  for all query spectra (eventually returning duplicated elements for
  query spectra matching multiple target spectra) and the values for the
  target spectra matched to the respective query spectra. See help on
  `$` above or examples below for details.

- `spectraVariables` returns all available spectra variables in the
  *query* and *target* spectra. The prefix `"target_"` is used to label
  spectra variables of target spectra (e.g. the name of the spectra
  variable for the MS level of target spectra is called
  `"target_msLevel"`).

- `target` returns the *target* `Spectra`.

- `targetVariables` returns the `spectraVariables` of *target* (prefixed
  with `"target_"`).

- `whichTarget` returns an `integer` with the indices of the spectra in
  *target* that match at least on spectrum in *query*.

- `whichQuery` returns an `integer` with the indices of the spectra in
  *query* that match at least on spectrum in *target*.

## Data manipulation and plotting

- `addProcessing`: add a processing step to both the *query* and
  *target* `Spectra` in `object`. Additional parameters for `FUN` can be
  passed *via* `...`. See `addProcessing` documentation in
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for more information.

- `plotSpectraMirror`: creates a mirror plot between the query and each
  matching target spectrum. Can only be applied to a `MatchedSpectra`
  with a single query spectrum. Setting parameter `scalePeaks = TRUE`
  will scale the peak intensities per spectrum to a total sum of one for
  a better graphical visualization. Additional plotting parameters can
  be passed through `...`. The parameters `ppm` and `tolerance` can be
  used to define the m/z tolerance for matching peaks between the query
  and target spectra. If not provided by the user, the values from the
  `param` object used to create the `MatchedSpectra` object are used; if
  these are missing, the default values (`ppm =20` and `tolerance = 0`)
  are used.

- `setBackend`: allows to change the *backend* of both the query and
  target
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object. The function will return a `MatchedSpectra` object with the
  query and target `Spectra` changed to the specified `backend`, which
  can be any backend extending
  [Spectra::MsBackend](https://rdrr.io/pkg/Spectra/man/MsBackend.html).

## See also

[`Matched()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/Matched.md)
for additional functions available for `MatchedSpectra`.

## Author

Johannes Rainer

## Examples

``` r

## Creating a dummy MatchedSpectra object.
library(Spectra)
df1 <- DataFrame(
    msLevel = 2L, rtime = 1:10,
    spectrum_id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"))
df2 <- DataFrame(
    msLevel = 2L, rtime = rep(1:10, 20),
    spectrum_id = rep(c("A", "B", "C", "D", "E"), 20))
sp1 <- Spectra(df1)
sp2 <- Spectra(df2)
## Define matches between query spectrum 1 with target spectra 2 and 5,
## query spectrum 2 with target spectrum 2 and query spectrum 4 with target
## spectra 8, 12 and 15.
ms <- MatchedSpectra(
    sp1, sp2, matches = data.frame(query_idx = c(1L, 1L, 2L, 4L, 4L, 4L),
                                   target_idx = c(2L, 5L, 2L, 8L, 12L, 15L),
                                   score = 1:6))

## Which of the query spectra match at least one target spectrum?
whichQuery(ms)
#> [1] 1 2 4

## Extracting spectra variables: accessor methods for spectra variables act
## as "left joins", i.e. they return a value for each query spectrum, with
## eventually duplicated elements if one query spectrum matches more than
## one target spectrum.

## Which target spectrum matches at least one query spectrum?
whichTarget(ms)
#> [1]  2  5  8 12 15

## Extracting the retention times of the query spectra.
ms$rtime
#>  [1]  1  1  2  3  4  4  4  5  6  7  8  9 10

## We have duplicated retention times for query spectrum 1 (matches 2 target
## spectra) and 4 (matches 3 target spectra). The retention time is returned
## for each query spectrum.

## Extracting retention times of the target spectra. Note that only retention
## times for target spectra matching at least one query spectrum are returned
## and an NA is reported for query spectra without matching target spectrum.
ms$target_rtime
#>  [1]  2  5  2 NA  8  2  5 NA NA NA NA NA NA

## The first query spectrum matches target spectra 2 and 5, thus their
## retention times are returned as well as the retention time of the second
## target spectrum that matches also query spectrum 2. The 3rd query spectrum
## does match any target spectrum, thus `NA` is returned. Query spectrum 4
## matches target spectra 8, 12, and 15, thus the next reported retention
## times are those from these 3 target spectra. None of the remaining 6 query
## spectra matches any target spectra and thus `NA` is reported for each of
## them.

## With `queryIndex` and `targetIndex` it is possible to extract the indices
## of the matched query-index pairs
queryIndex(ms)
#> [1] 1 1 2 4 4 4
targetIndex(ms)
#> [1]  2  5  2  8 12 15

## The first match is between query index 1 and target index 2, the second
## match between query index 1 and target index 5 and so on.
## We could use these indices to extract a `Spectra` object containing only
## matched target spectra and assign a spectra variable with the indices of
## the query spectra
matched_target <- target(ms)[targetIndex(ms)]
matched_target$query_index <- queryIndex(ms)

## This `Spectra` object thus contains information from the matching, but
## is a *conventional* `Spectra` object that could be used for further
## analyses.

## `spectraData` can be used to extract all (or selected) spectra variables
## from the object. Same as with `$`, a left join between the specta
## variables from the query spectra and the target spectra is performed. The
## prefix `"target_"` is used to label the spectra variables from the target
## spectra. Below we extract selected spectra variables from the object.
res <- spectraData(ms, columns = c("rtime", "spectrum_id",
    "target_rtime", "target_spectrum_id"))
res
#> DataFrame with 13 rows and 4 columns
#>         rtime spectrum_id target_rtime target_spectrum_id
#>     <integer> <character>    <integer>        <character>
#> 1           1           a            2                  B
#> 2           1           a            5                  E
#> 3           2           b            2                  B
#> 4           3           c           NA                 NA
#> 5           4           d            8                  C
#> ...       ...         ...          ...                ...
#> 9           6           f           NA                 NA
#> 10          7           g           NA                 NA
#> 11          8           h           NA                 NA
#> 12          9           i           NA                 NA
#> 13         10           j           NA                 NA
res$spectrum_id
#>  [1] "a" "a" "b" "c" "d" "d" "d" "e" "f" "g" "h" "i" "j"
res$target_spectrum_id
#>  [1] "B" "E" "B" NA  "C" "B" "E" NA  NA  NA  NA  NA  NA 

## Again, all values for query spectra are returned and for query spectra not
## matching any target spectrum NA is reported as value for the respecive
## variable.

## The example matched spectra object contains all query and all target
## spectra. Below we subset the object keeping only query spectra that are
## matched to at least one target spectrum.
ms_sub <- ms[whichQuery(ms)]

## ms_sub contains now only 3 query spectra:
length(query(ms_sub))
#> [1] 3

## while the original object contains all 10 query spectra:
length(query(ms))
#> [1] 10

## Both object contain however still the full target `Spectra`:
length(target(ms))
#> [1] 200
length(target(ms_sub))
#> [1] 200

## With the `pruneTarget` we can however reduce also the target spectra to
## only those that match at least one query spectrum
ms_sub <- pruneTarget(ms_sub)
length(target(ms_sub))
#> [1] 5
```
