# MetaboAnnotation 1.5

## Changes in 1.5.7

- Add function `.group_standards_iteration` to allow iteration through matrix 
  of standards and group them if they are dissimilar enough.

## Changes in 1.5.6

- Fix issue in the vignette. Thanks @RemyDeB for the fix.

## Changes in 1.5.5

- Update objects to the new definitions in `Spectra` version 1.11.10.

## Changes in 1.5.4

- Add functions `targetIndex` and `queryIndex` to extract the indices of the
  matched pairs query-target.
- Add examples and a section to the vignette explaining their use.

## Changes in 1.5.3

- Add support to `matchValues` for matching between `data.frame` and `Spectra`
  objects.

## Changes in 1.5.2

- Fix vignette, examples and unit tests using `QFeatures`.
- Import `query` from `AnnotationHub`.

## Changes in 1.5.1

- Add possibility to select the spectra variable for retention time matching in
  `matchSpectra` (issue
  [#98](https://github.com/rformassspectrometry/MetaboAnnotation/issues/98)).

# MetaboAnnotation 1.3

## Changes in 1.3.2

- Add `mzR` as suggested package to ensure package vignettes can be built.

## Changes in 1.3.1

- Small changes in `matchSpectra` to avoid unnecessary object creation.
- Use `backendBpparam` to disable parallel processing of `matchSpectra` if the
  backend does not support it.

# MetaboAnnotation 1.1

## Changes in 1.1.6

- `scoreVariables` function to return the names of the score variables in a
  `Matched` object.

## Changes in 1.1.5

- Fix issues on BioC build machines.

## Changes in 1.1.4

- `matchSpectra`: support a `CompDb` with parameter `target`.
- Add `CompAnnotionSource` classes to support definition of references to
  annotation resources.
- Add `CompDbSource` class defining a reference to a `CompDb` database.
- `matchSpectra`: support for `CompDbSource` with parameter `target`.

## Changes in 1.1.3

- Extend `filterMatches` framework (issue
  [#86](https://github.com/rformassspectrometry/MetaboAnnotation/issues/86)).
  `ScoreThresholdParam` added to perform filtering the matches based on a
  threshold for the scores.

## Changes in 1.1.2

- `lapply` and `endoapply` methods (issue
  [#84](https://github.com/rformassspectrometry/MetaboAnnotation/issues/84)).
  `lapply` allows to apply any function to each subset of matches for each
  `query` element and returns the corresponding `list` of results.
  `endoapply` is similar but applies a function returning a `Matched` and
  returns a `Matched` representing updated matches.

## Changes in 1.1.1

- Extend `filterMatches` framework (issue
  [#81](https://github.com/rformassspectrometry/MetaboAnnotation/issues/81)).
  `SelectMatchesParam` and `TopRankedMatchesParam` added to perform respectively
  manual filtering and keeping only the best ranked matches for each `query`
  element.

# MetaboAnnotation 0.99

## Changes in 0.99.15

- Highlight query and target spectra in different colors for
  `validateMatchedSpectra`.
- `query` and/or `target` of type `SummarizedExperiment` supported for
  `Matched` objects.
- `MatchedSummarizedExperiment` class removed.
- `query` and/or `target` of type `QFeatures` supported for `Matched` objects.
- Support `SummarizedExperiment` and `QFeatures` for both `query` and `target`
  parameters in `matchValues`.

## Changes in 0.99.14

- Improve plotly-based mirror plots in `validateMatchedSpectra`.

## Changes in 0.99.13

- Fix issue about `matchedData` not working for result objects of
  `matchValues, Mz2MassParam` and `matchValues, Mz2MassRtParam` (issue
  [#69](https://github.com/rformassspectrometry/MetaboAnnotation/issues/69)).

## Changes in 0.99.12

- Update plotly-based mirror plots in `validateMatchedSpectra`.

## Changes in 0.99.11

- Change `matchMz` into `matchValues` (issue
  [#65](https://github.com/rformassspectrometry/MetaboAnnotation/issues/65)).

## Changes in 0.99.10

- Add `validateMatchedSpectra` for manual inspection and validation of an
  `MatchedSpectra` object.

## Changes in 0.99.9

- Add `setBackend` for `MatchedSpectra` objects.

## Changes in 0.99.8

- Add `matchMz, Mz2MassParam` and `matchMz, Mz2MassRtParam`. (issue
  [#56](https://github.com/rformassspectrometry/MetaboAnnotation/issues/56)).

## Changes in 0.99.7

- Add formula matching functions.

## Changes in 0.99.5

- Add parameter `...` to `plotSpectraMirror`.
- Definitions of "`score`", "`score_rt`" changed to be the difference
  (with sign) between query and target m/z or retention time respectively.
- `"ppm_error"` becomes error without sign.

## Changes in 0.99.4

- Add matches m/z error (variable `"ppm_error"`) to the `Matched` object
  returned by `matchMz`.

## Changes in 0.99.3

- Address Herve's comments.


# MetaboAnnotation 0.2

## Changes in 0.2.11

- Fix calculation of correct number of rows/columns of the plot in
  `plotSpectraMirror`.

## Changes in 0.2.10

- Add parameter `toleranceRt` to `CompareSpectraParam` to enable retention
  time-based pre-filtering (issue
  [#35](https://github.com/rformassspectrometry/MetaboAnnotation/issues/35)).

## Changes in 0.2.9

- Add support for manually defined adducts to `Mass2MzParam` (issue
  [#41](https://github.com/rformassspectrometry/MetaboAnnotation/issues/41)).

## Changes in 0.2.8

- Add parameter `THRESHFUN_REVERSE` to `MatchForwardReverseParam` to allow
  filtering results on forward **and** reverse score (issue
  [#37](https://github.com/rformassspectrometry/MetaboAnnotation/issues/37)).

## Changes in 0.2.7

- Performance improvement in `matchSpectra` if no precursor m/z filter is used
  (issue
  [#38](https://github.com/rformassspectrometry/MetaboAnnotation/issues/38)).
- Report number of matching peaks in `matchSpectra,MatchForwardReverseParam`
  (issue
  [#36](https://github.com/rformassspectrometry/MetaboAnnotation/issues/36)).

## Changes in 0.2.6

- Fix bug in `matchSpectra` that was wrongly calculating the acceptable m/z
  difference if `tolerance` was > 0 (issue
  [#34](https://github.com/rformassspectrometry/MetaboAnnotation/issues/34)).
  Fix proposed by Hugo Varet (@hvaret).

## Changes in 0.2.5

- Improve performance of `matchMz`.
- Rename `queryColumn` and `targetColumn` to `queryColname` and `targetColname`.

## Changes in 0.2.4

- Support `data.frame`, `DataFrame` and `matrix` in `matchMz`.
- Add `addMatches` and `filterMatches` functions.

## Changes in 0.2.3

- Fixes in `MatchedSpectra`.

## Changes in 0.2.2

- Add `MatchedSummarizedExperiment`.

## Changes in 0.2.1

- Rename `TargetMass2MzParam` to `Mass2MzParam`.

## Changes in 0.2.0

- Add support for matching m/z against m/z and m/z in addition to retention
  times to `matchMz`.

# MetaboAnnotation 0.0

## Changes in 0.0.4

- Fix vignette, documentations and unit tests.
