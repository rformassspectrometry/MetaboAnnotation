# Spectral matching

The `matchSpectra` method matches (compares) spectra from `query` with
those from `target` based on settings specified with `param` and returns
the result from this as a
[MatchedSpectra](https://rformassspectrometry.github.io/MetaboAnnotation/reference/MatchedSpectra.md)
object.

## Usage

``` r
matchSpectra(query, target, param, ...)
```

## Arguments

- query:

  [Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object with the (experimental) spectra.

- target:

  MS data to compare against. Can be another
  [Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html).

- param:

  parameter object containing the settings for the matching (e.g.
  eventual prefiltering settings, cut-off value for similarity above
  which spectra are considered matching etc).

- ...:

  optional parameters.

## Value

a
[MatchedSpectra](https://rformassspectrometry.github.io/MetaboAnnotation/reference/MatchedSpectra.md)
object with the spectra matching results.

## See also

[`CompareSpectraParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompareSpectraParam.md)
for the comparison between
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html)
objects.

## Author

Johannes Rainer
