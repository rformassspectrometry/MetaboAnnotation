#' @include CompAnnotationSource.R

#' @title Compond Annotation Sources for databases supporting `Spectra`
#'
#' @aliases SpectraDbSource-class
#'
#' @description
#'
#' `SpectraDbSource` objects represent references to annotation databases that
#' are supported by a dedicated [MsBackend] backend and hence allow to return
#' their data as a [Spectra] object. A specific `MsBackend` has thus to be
#' available (and needs to be supported) by the database. Examples are
#' *MassBank* with the
#' [MsBackendMassbank](https://github.com/RforMassSpectrometry/MsBackendMassbank)
#' backend or *WeizMass* with the
#' [MsBackendWeizMass](https://github.com/RforMassSpectrometry/MsBackendWeizMass).
#' Additional parameters that need to be defined are the type of SQL database
#' system that is used by the resource (e.g. `MariaDB` for MySQL/MariaDB
#' servers) as well as credential information allowing to connect to the
#' database. Connection to the database is only performed in the `matchSpectra`
#' call and the returned [MatchedSpectra] object contains only matching spectra
#' from the reference (target) database.
#'
#' New `SpectraDbSource` objects should not be created manually, but using one
#' of the constructor functions listed below. Functions are available for
#' different annotation databases that provide also useful defaults simplifying
#' thus access to these resources.
#'
#' - `RemoteWeizMassSource`: connects to a *remote* WeizMass (MySQL) database
#'   hosted on a remote server. Parameter `user`, `pass`, `host` and `dbname`
#'   are optional and allow connection to a non-standard (default) WeizMass
#'   database (e.g. a local installation).
#'   For WeizMass databases, users have to agree to the WeizMass license by
#'   specifically setting parameter `IagreeToTheLicense` to `TRUE`.
#'
#' - `LocalWeizMassSource`: connects to a *local* WeizMass database available as
#'   a *SQLite* database file. The database file has to be provided with the
#'   `dbfile` parameter.
#'   For WeizMass databases, users have to agree to the WeizMass license by
#'   specifically setting parameter `IagreeToTheLicense` to `TRUE`.
#'
#' @param dbfile For `LocalWeizMassSource`: `character(1)` with the (SQLite)
#'     database file.
#'
#' @param dbname For `RemoteWeizMassSource`: optional `character(1)` with the
#'     name of the database on the server.
#'
#' @param host For `RemoteWeizMassSource`: optional `character(1)` with the
#'     name (address) of the host with the WeizMass database.
#'
#' @param IagreeToTheLicense `logical(1)` specifying whether the user accepts
#'     to the WeizMass license agreement or not.
#'
#' @param object A `SpectraDbSource` object.
#'
#' @param pass For `RemoteWeizMassSource`: optional `character(1)` with the
#'     password to access the WeizMass database.
#'
#' @param user For `RemoteWeizMassSource`: optional `character(1)` with the
#'     username to access the WeizMass database.
#'
#' @author Johannes Rainer, Nir Shachaf
#'
#' @family Compound annotation sources
#'
#' @name SpectraDbSource
#'
NULL

#' @importClassesFrom DBI DBIDriver
#'
#' @importClassesFrom Spectra MsBackend
setClass(
    "SpectraDbSource",
    contains = "CompAnnotationSource",
    slots = c(host = "character",
              user = "character",
              pass = "character",
              dbname = "character",
              description = "character",
              drv = "DBIDriver",
              backend = "MsBackend"),
    prototype = list(host = character(),
                     user = character(),
                     pass = character(),
                     dbname = character(),
                     description = character(),
                     drv = NULL,
                     backend = NULL))

#' @export
#'
#' @rdname SpectraDbSource
RemoteWeizMassSource <- function(IagreeToTheLicense = FALSE, host = character(),
                                 user = character(), pass = character(),
                                 dbname = character()) {
    if (!IagreeToTheLicense)
        stop("You must agree to the WeizMass license agreement.")
    if (!requireNamespace("MsBackendWeizMass", quietly = TRUE))
        stop("The use of 'RemoteWeizMassSource' requires package ",
             "'MsBackendWeizMass'. Please install it with ",
             "'BiocInstaller::install(",
             "\"RforMassSpectrometry/MsBackendWeizMass\")'.")
    if (!requireNamespace("RMariaDB", quietly = TRUE))
        stop("The use of 'RemoteWeizMassSource' requires package ",
             "'RMariaDB'. Please install it with ",
             "'BiocInstaller::install(\"RMariaDB\")'.")
    backend <- MsBackendWeizMass::MsBackendWeizMass()
    drv <- RMariaDB::MariaDB()
    new("SpectraDbSource", host = host, user = user, pass = pass,
        dbname = dbname, drv = drv, backend = backend,
        description = "Spectra source: WeizMass")
}

#' @export
#'
#' @rdname SpectraDbSource
LocalWeizMassSource <- function(IagreeToTheLicense = FALSE,
                                dbfile = character()) {
    if (!IagreeToTheLicense)
        stop("You must agree to the WeizMass license agreement.")
    if (!requireNamespace("MsBackendWeizMass", quietly = TRUE))
        stop("The use of 'LocalWeizMassSource' requires package ",
             "'MsBackendWeizMass'. Please install it with ",
             "'BiocInstaller::install(",
             "\"RforMassSpectrometry/MsBackendWeizMass\")'.")
    if (!requireNamespace("RSQLite", quietly = TRUE))
        stop("The use of 'LocalWeizMassSource' requires package ",
             "'RSQLite'. Please install it with ",
             "'BiocInstaller::install(\"RSQLite\")'.")
    backend <- MsBackendWeizMass::MsBackendWeizMass()
    drv <- RSQLite::SQLite()
    new("SpectraDbSource", dbname = dbfile, drv = drv, backend = backend,
        description = paste0("Spectra source: WeizMass\n",
                             "Database file: ", basename(dbfile)))
}

## #' @export
## #'
## #' @rdname SpectraDbSource
## MassBankDbSource <- function(dbname = "MassBank", host = "localhost",
##                              user = "massbank", pass = "massbank",
##                              drv = RMariaDB::MariaDB()) {
##     if (!requireNamespace("MsBackendMassbank", quietly = TRUE))
##         stop("The use of 'MassBankDbSource' requires package ",
##              "'MsBackendMassbank'. Please install it with ",
##              "'BiocInstaller::install(\"MsBackendMassbank\")'.")
##     new("SpectraDbSource", dbname = dbname, drv = drv,
##         backend = MsBackendMassbank::MsBackendMassbankSql(),
##         user = user, host = host, pass = pass)
## }

#' @rdname SpectraDbSource
setMethod("show", "SpectraDbSource", function(object) {
    callNextMethod()
    cat(object@description, "\n")
    cat("MsBackend:", class(object@backend)[[1L]], "\n")
})

#' @importFrom utils read.table
.get_weizmass_conf <- function() {
    ## Options:
    ## - read from file within the package
    ## - get from environment variables
    ## - define within the source code
    cfg <- read.table(system.file("cfg", "wm.cfg",
                                  package = "MetaboAnnotation"),
                      sep = "=")
    c(user = cfg[cfg[, 1] == "user", 2],
      pass = cfg[cfg[, 1] == "pass", 2],
      host = cfg[cfg[, 1] == "host", 2],
      dbname = cfg[cfg[, 1] == "dbname", 2])
}

#' @rdname CompareSpectraParam
#'
#' @importFrom Spectra MsBackendDataFrame
#'
#' @importMethodsFrom Spectra backendInitialize
#'
#' @importClassesFrom Spectra MsBackendDataFrame
#'
#' @importMethodsFrom DBI dbConnect
setMethod(
    "matchSpectra", signature(query = "Spectra", target = "SpectraDbSource",
                              param = "Param"),
    function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        if (inherits(target@drv, "SQLiteDriver")) {
            con <- dbConnect(target@drv, target@dbname)
        } else {
            if (!length(target@user)) {
                if (inherits(target@backend, "MsBackendWeizMass")) {
                    cr <- .get_weizmass_conf()
                    target@user <- cr["user"]
                    target@pass <- cr["pass"]
                } else stop("Missing 'user' and 'pass'")
            }
            con <- dbConnect(target@drv, host = target@host,
                             user = target@user, pass = target@pass,
                             dbname = target@dbname)
        }
        ## create a specific backend
        be <- backendInitialize(target@backend, con)
        ## get the Spectra from the backend and call matchSpectra
        res <- matchSpectra(query, Spectra(be), param = param,
                            BPPARAM = BPPARAM)
        ## keep only matching reference/target spectra and change the
        ## backend to MsBackendDataFrame
        res <- pruneTarget(res)
        res@target <- setBackend(res@target, backend = MsBackendDataFrame())
        res
    })
