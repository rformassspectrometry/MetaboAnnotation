#' @include CompAnnotationSource.R

#' @title Compound Annotation Sources for `CompDb` databases
#'
#' @aliases CompDbSource-class
#'
#' @description
#'
#' `CompDbSource` objects represent references to [CompDb] database-backed
#' annotation resources. Instances are expected to be created with the dedicated
#' construction functions such as `MassBankSource` or the generic
#' `CompDbSource`. The annotation data is not stored within the object but will
#' be accessed/loaded within the object's `matchSpectra` method.
#'
#' New `CompDbSource` objects can be created using the functions:
#'
#' - `CompDbSource`: create a new `CompDbSource` object from an existing
#'   `CompDb` database. The (SQLite) database file (including the full path)
#'   needs to be provided with parameter `dbfile`.
#'
#' - `MassBankSource`: retrieves a `CompDb` database for the specified MassBank
#'   release from Bioconductor's online `AnnotationHub` (if it exists) and
#'   uses that. Note that `AnnotationHub` resources are cached locally and thus
#'   only downloaded the first time.
#'   The function has parameters `release` which allows to define the desired
#'   MassBank release (e.g. `release = "2021.03"` or `release = "2022.06"`)
#'   and `...` which allows to pass optional parameters to the `AnnotationHub`
#'   constructor function, such as `localHub = TRUE` to use only the cached
#'   data and avoid updating/retrieving updates from the internet.
#'
#' Other functions:
#'
#' - `metadata`: get metadata (information) on the annotation resource.
#'
#' @param dbfile `character(1)` with the database file (including the full
#'     path).
#'
#' @param object A `CompDbSource` object.
#'
#' @param release A `character(1)` defining the version/release of MassBank that
#'     should be used.
#'
#' @param x A `CompDbSource` object.
#'
#' @param ... For `CompDbSource`: ignored. For `MassBankSource`: optional
#'     parameters passed to the `AnnotationHub` constructor function.
#'
#' @name CompDbSource
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Locate a CompDb SQLite database file. For this example we use the test
#' ## database from the `CompoundDb` package.
#' fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
#' ann_src <- CompDbSource(fl)
#'
#' ## The object contains only the reference/link to the annotation resource.
#' ann_src
#'
#' ## Retrieve a CompDb with MassBank data for a certain MassBank release
#' mb_src <- MassBankSource("2021.03")
#' mb_src
NULL

#' @description
#'
#' It is better to NOT put neither a connection object nor the data itself into
#' the source object to allow parallel processing. The `matchSpectra` method
#' needs then to make sure to load/retrieve the data and connect to it.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @exportClass CompDbSource
setClass(
    "CompDbSource",
    contains = "CompAnnotationSource",
    slots = c(dbfile = "character"),
    prototype = list(dbfile = character())
)

#' validator to check if the provided file name is indeed a CompDb database.
#'
#' @importFrom CompoundDb CompDb
#'
#' @noRd
.validate_dbfile <- function(x) {
    if (length(x))
        validObject(CompDb(x))
    else TRUE
}

setValidity("CompDbSource", function(object) {
    .validate_dbfile(object@dbfile)
})

#' @rdname CompDbSource
#'
#' @export
CompDbSource <- function(dbfile = character()) {
    new("CompDbSource", dbfile = dbfile)
}

#' @export
#'
#' @rdname CompDbSource
setMethod("metadata", "CompDbSource", function(x, ...) {
    db <- CompDb(x@dbfile)
    metadata(db)
})

#' @importFrom methods callNextMethod
#'
#' @rdname CompDbSource
setMethod("show", "CompDbSource", function(object) {
    callNextMethod()
    md <- metadata(object)
    cat("Metadata information:\n")
    cat(paste0("  - ", md[, 1], ": ", md[, 2], "\n"), sep = "")
})

#' @rdname CompareSpectraParam
#'
#' @importFrom Spectra MsBackendDataFrame
#'
#' @importClassesFrom Spectra MsBackendDataFrame
setMethod(
    "matchSpectra", signature(query = "Spectra", target = "CompDbSource",
                              param = "Param"),
    function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        ## connect to the database
        db <- CompDb(target@dbfile)
        ## get the Spectra from the source and call matchSpectra
        res <- matchSpectra(query, Spectra(db), param = param,
                            BPPARAM = BPPARAM)
        ## keep only matching reference/target spectra  and change the
        ## backend to MsBackendDataFrame
        res <- pruneTarget(res)
        res@target <- setBackend(res@target, backend = MsBackendDataFrame())
        res
    })

#' @export
#'
#' @rdname CompDbSource
MassBankSource <- function(release = "2021.03", ...) {
    if (!requireNamespace("AnnotationHub", quietly = TRUE))
        stop("'MassBankSource' requires the 'AnnotationHub' package.\n",
             "Please install it with 'BiocManager::install(",
             "\"AnnotationHub\")' and try again.")
    ah <- AnnotationHub::AnnotationHub(...)
    res <- AnnotationHub::query(ah, c("MassBank", release))
    if (!length(res))
        stop("MassBank release \"", release, "\" not found in AnnotationHub")
    if (length(res) > 1)
        stop("Provided release is ambiguous: ", length(res), " data sets in ",
             "AnnotationHub match the provided release information.")
    fn <- unname(AnnotationHub::cache(res))
    CompDbSource(fn)
}
