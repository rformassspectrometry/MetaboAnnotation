# Validating MatchedSpectra

The `validateMatchedSpectra()` function opens a simple shiny application
that allows to browse results stored in a `MatchedSpectra` object and to
*validate* the presented matches. For each query spectrum a table with
matched target spectra are shown (if available) and an interactive
mirror plot is generated. Valid matches can be selected using a check
box which is displayed below the mirror plot. Upon pushing the "Save &
Close" button the app is closed and a filtered `MatchedSpectra` is
returned, containing only *validated* matches.

Note that column `"query_index_"` and `"target_index_"` are temporarily
added to the query and target `Spectra` object to display them in the
interactive graphics for easier identification of the compared spectra.

## Usage

``` r
validateMatchedSpectra(object)
```

## Arguments

- object:

  A non-empty instance of class `MatchedSpectra`.

## Value

A `MatchedSpectra` with validated results.

## Author

Carolin Huber, Michael Witting, Johannes Rainer

## Examples

``` r

library(Spectra)
## Load test data from *MsDataHub*
fl <- MsDataHub::PestMix1_DDA.mzML()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))

## Normalize intensities and match spectra
csp <- CompareSpectraParam(requirePrecursor = TRUE,
                           THRESHFUN = function(x) x >= 0.7)
norm_int <- function(x) {
    x[, "intensity"] <- x[, "intensity"] / max(x[, "intensity"]) * 100
    x
}
ms <- matchSpectra(addProcessing(pest_ms2, norm_int),
                   addProcessing(minimb, norm_int), csp)

## validate matches using the shiny app. Note: the call is only executed
## in interactive mode.
if (interactive()) {
    res <- validateMatchedSpectra(ms)
}
```
