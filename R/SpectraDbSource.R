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
#' access to these resources.
#'
#' - `WeizMassSource`: provides access to a WeizMass database. This database
#'   can be a SQLite database (in which case the `sqlite` parameter has to be
#'   set to `TRUE` and the database file has to be provided with the `dbname`
#'   parameter) or a (remote) WeizMass database. In the latter case the version
#'   of the resource has to be provided. Alternatively, using parameters `host`,
#'   `dbname`, `user` and `pass` it is possible to connect to a custom location
#'   of a WeizMass database.
#'   For WeizMass databases, users have to agree to the WeizMass license by
#'   specifically setting parameter `IagreeToTheLicense` to `TRUE`.
#'
#' @param dbname For `WeizMassSource`: optional `character(1)` with the
#'     name of the database on the server. If `sqlite = TRUE` the file name
#'     of the SQLite database needs to be defined with this parameter.
#'
#' @param host For `WeizMassSource`: optional `character(1)` with the
#'     name (address) of the host with the WeizMass database.
#'
#' @param IagreeToTheLicense `logical(1)` specifying whether the user accepts
#'     to the WeizMass license agreement or not.
#'
#' @param object A `SpectraDbSource` object.
#'
#' @param pass For `WeizMassSource`: optional `character(1)` with the
#'     password to access the WeizMass database.
#'
#' @param user For `WeizMassSource`: optional `character(1)` with the
#'     username to access the WeizMass database.
#'
#' @param version For `WeizMassSource`: `character(1)` specifying the version
#'     (release) of the database to connect to.
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
              version = "character",
              description = "character",
              drv = "DBIDriver",
              backend = "MsBackend"),
    prototype = list(host = character(),
                     user = character(),
                     pass = character(),
                     dbname = character(),
                     version = character(),
                     description = character(),
                     drv = NULL,
                     backend = NULL))

.sql_weiz_mass <- function(host = character(), user = character(),
                           pass = character(), dbname = character(),
                           version = character()) {
    if (!requireNamespace("RMariaDB", quietly = TRUE))
        stop("The use of 'RemoteWeizMassSource' requires package ",
             "'RMariaDB'. Please install it with ",
             "'BiocInstaller::install(\"RMariaDB\")'.")
    backend <- MsBackendWeizMass::MsBackendWeizMass()
    drv <- RMariaDB::MariaDB()
    vrs <- ""
    if (length(version))
        vrs <- paste0(" version ", version)
    new("SpectraDbSource", host = host, user = user, pass = pass,
        dbname = dbname, drv = drv, backend = backend, version = version,
        description = paste0("Spectra source: WeizMass", vrs))
}

.sqlite_weiz_mass <- function(dbname = character(), ...) {
    if (!requireNamespace("RSQLite", quietly = TRUE))
        stop("The use of 'LocalWeizMassSource' requires package ",
             "'RSQLite'. Please install it with ",
             "'BiocInstaller::install(\"RSQLite\")'.")
    if (!length(dbname))
        stop("Parameter 'dbname' is mandatory for 'sqlite = TRUE'")
    backend <- MsBackendWeizMass::MsBackendWeizMass()
    drv <- RSQLite::SQLite()
    new("SpectraDbSource", dbname = dbname, drv = drv, backend = backend,
        description = paste0("Spectra source: WeizMass\n",
                             "Database file: ", basename(dbname)))
}

#' @export
#'
#' @rdname SpectraDbSource
WeizMassSource <- function(IagreeToTheLicense = FALSE, sqlite = FALSE,
                           version = character(), user = character(),
                           pass = character(), host = character(),
                           dbname = character()) {
    if (!IagreeToTheLicense)
        stop("You must agree to the WeizMass license agreement.")
    if (!requireNamespace("MsBackendWeizMass", quietly = TRUE))
        stop("The use of 'LocalWeizMassSource' requires package ",
             "'MsBackendWeizMass'. Please install it with ",
             "'BiocInstaller::install(",
             "\"RforMassSpectrometry/MsBackendWeizMass\")'.")
    if (sqlite)
        .sqlite_weiz_mass(dbname = dbname)
    else .sql_weiz_mass(version = version, user = user, pass = pass,
                        host = host, dbname = dbname)
}

#' @rdname SpectraDbSource
setMethod("show", "SpectraDbSource", function(object) {
    callNextMethod()
    cat(object@description, "\n")
    cat("MsBackend:", class(object@backend)[[1L]], "\n")
})

#' @importFrom utils read.table
.get_weizmass_conf <- function(version = character()) {
    ## Options:
    ## - read from file within the package
    ## - get from environment variables
    ## - define within the source code
    fl <- system.file("cfg", paste0("wm", version, ".cfg"),
                      package = "MetaboAnnotation")
    if (fl == "")
        stop("No configuration for WeizMass version \"", version, "\" found")
    cfg <- read.table(fl, sep = "=")
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
                    cr <- .get_weizmass_conf(target@version)
                    target@user <- cr["user"]
                    target@pass <- cr["pass"]
                    target@host <- cr["host"]
                    target@dbname <- cr["dbname"]
                } else stop("Missing parameters 'dbname', 'host', 'user'",
                            " and 'pass'")
            }
            con <- dbConnect(target@drv, host = target@host,
                             user = target@user, pass = target@pass,
                             dbname = target@dbname)
        }
        be <- backendInitialize(target@backend, con)
        res <- matchSpectra(query, Spectra(be), param = param,
                            BPPARAM = BPPARAM)
        ## keep only matching reference/target spectra and change the
        ## backend to MsBackendDataFrame
        res <- pruneTarget(res)
        res@target <- setBackend(res@target, backend = MsBackendDataFrame())
        res
    })
