# Compound Annotation Sources for `CompDb` databases

`CompDbSource` objects represent references to
[CompoundDb::CompDb](https://rdrr.io/pkg/CompoundDb/man/CompDb.html)
database-backed annotation resources. Instances are expected to be
created with the dedicated construction functions such as
`MassBankSource` or the generic `CompDbSource`. The annotation data is
not stored within the object but will be accessed/loaded within the
object's `matchSpectra` method.

New `CompDbSource` objects can be created using the functions:

- `CompDbSource`: create a new `CompDbSource` object from an existing
  `CompDb` database. The (SQLite) database file (including the full
  path) needs to be provided with parameter `dbfile`.

- `MassBankSource`: retrieves a `CompDb` database for the specified
  MassBank release from Bioconductor's online `AnnotationHub` (if it
  exists) and uses that. Note that `AnnotationHub` resources are cached
  locally and thus only downloaded the first time. The function has
  parameters `release` which allows to define the desired MassBank
  release (e.g. `release = "2021.03"` or `release = "2022.06"`) and
  `...` which allows to pass optional parameters to the `AnnotationHub`
  constructor function, such as `localHub = TRUE` to use only the cached
  data and avoid updating/retrieving updates from the internet.

Other functions:

- `metadata`: get metadata (information) on the annotation resource.

## Usage

``` r
CompDbSource(dbfile = character())

# S4 method for class 'CompDbSource'
metadata(x, ...)

# S4 method for class 'CompDbSource'
show(object)

MassBankSource(release = "2021.03", ...)
```

## Arguments

- dbfile:

  `character(1)` with the database file (including the full path).

- x:

  A `CompDbSource` object.

- ...:

  For `CompDbSource`: ignored. For `MassBankSource`: optional parameters
  passed to the `AnnotationHub` constructor function.

- object:

  A `CompDbSource` object.

- release:

  A `character(1)` defining the version/release of MassBank that should
  be used.

## Author

Johannes Rainer

## Examples

``` r

## Locate a CompDb SQLite database file. For this example we use the test
## database from the `CompoundDb` package.
fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
ann_src <- CompDbSource(fl)

## The object contains only the reference/link to the annotation resource.
ann_src
#> Object of class CompDbSource 
#> Metadata information:
#>   - source: MassBank
#>   - url: https://massbank.eu/MassBank/
#>   - source_version: 2020.09
#>   - source_date: 1603272565
#>   - organism: NA
#>   - db_creation_date: Thu Oct 22 08:45:31 2020
#>   - supporting_package: CompoundDb
#>   - supporting_object: CompDb

## Retrieve a CompDb with MassBank data for a certain MassBank release
mb_src <- MassBankSource("2021.03")
mb_src
#> Object of class CompDbSource 
#> Metadata information:
#>   - source: MassBank
#>   - url: https://massbank.eu/MassBank/
#>   - source_version: 2021.03
#>   - source_date: 2021-02-26
#>   - organism: NA
#>   - db_creation_date: Tue Aug 30 06:53:53 2022
#>   - supporting_package: CompoundDb
#>   - supporting_object: CompDb
```
