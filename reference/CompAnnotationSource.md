# Compound Annotation Sources

`CompAnnotationSource`s (i.e. classes extending the base virtual
`CompAnnotationSource` class) define and provide access to a
(potentially remote) compound annotation resource. This aims to simplify
the integration of external annotation resources by automating the
actual connection (or data resource download) process from the user. In
addition, since the reference resource is not directly exposed to the
user it allows integration of annotation resources that do not allow
access to the full data.

Objects extending `CompAnnotationSource` available in this package are:

- [`CompDbSource()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompDbSource.md):
  annotation source referencing an annotation source in the
  `[CompoundDb::CompDb()]` format ( from the `CompoundDb` Bioconductor
  package).

Classes extending `CompAnnotationSource` need to implement the
`matchSpectra` method with parameters `query`, `target` and `param`
where `query` is the `Spectra` object with the (experimental) query
spectra, `target` the object extending the `CompAnnotationSource` and
`param` the parameter object defining the similarity calculation (e.g.
[`CompareSpectraParam()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompareSpectraParam.md).
The method is expected to return a
[MatchedSpectra](https://rformassspectrometry.github.io/MetaboAnnotation/reference/MatchedSpectra.md)
object.

`CompAnnotationSource` objects are not expected to contain any
annotation data. Access to the annotation data (in form of a `Spectra`
object) is suggested to be only established within the object's
`matchSpectra` method. This would also enable parallel processing of
annotations as no e.g. database connection would have to be shared
across processes.

## Usage

``` r
# S4 method for class 'Spectra,CompAnnotationSource,Param'
matchSpectra(query, target, param, ...)

# S4 method for class 'CompAnnotationSource'
show(object)

# S4 method for class 'CompAnnotationSource'
metadata(x, ...)
```

## Arguments

- query:

  for `matchSpectra`:
  [Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object with the query spectra.

- target:

  for `matchSpectra`: object extending CompAnnotationSource (such as
  [CompDbSource](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompDbSource.md))
  with the target (reference) spectra to compare `query` against.

- param:

  for `matchSpectra`: parameter object (such as
  [CompareSpectraParam](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompareSpectraParam.md))
  defining the settings for the matching.

- ...:

  additional parameters passed to `matchSpectra`.

- object:

  A `CompAnnotationSource` object.

- x:

  A `CompAnnotationSource` object.

## Methods that need to be implemented

For an example implementation see
[`CompDbSource()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/CompDbSource.md).

- `matchSpectra`: function to match experimental MS2 spectra against the
  annotation source. See
  [`matchSpectra()`](https://rformassspectrometry.github.io/MetaboAnnotation/reference/matchSpectra.md)
  for parameters.

- `metadata`: function to provide metadata on the annotation resource
  (host, source, version etc).

- `show` (optional): method to provide general information on the data
  source.

## Author

Johannes Rainer, Nir Shachaf
