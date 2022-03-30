#' @title Validating MatchedSpectra
#'
#' @description
#'
#' The `validateMatchedSpectra()` function opens a simple shiny application
#' that allows to browse results stored in a `MatchedSpectra` object and to
#' *validate* the presented matches. For each query spectrum a table with
#' matched target spectra are shown (if available) and an interactive mirror
#' plot is generated. Valid matches can be selected using a check box which is
#' displayed below the mirror plot. Upon pushing the "Save & Close"
#' button the app is closed and a filtered `MatchedSpectra` is returned,
#' containing only *validated* matches.
#'
#' @param object A non-empty instance of class `MatchedSpectra`.
#'
#' @return A `MatchedSpectra` with validated results.
#'
#' @export
#'
#' @importFrom BiocParallel bpparam SerialParam register
#'
#' @author Carolin Huber, Michael Witting, Johannes Rainer
#'
#' @examples
#'
#' library(Spectra)
#' ## Load test data
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
#' pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
#' pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
#' load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))
#'
#' ## Normalize intensities and match spectra
#' csp <- CompareSpectraParam(requirePrecursor = TRUE,
#'                            THRESHFUN = function(x) x >= 0.7)
#' norm_int <- function(x) {
#'     x[, "intensity"] <- x[, "intensity"] / max(x[, "intensity"]) * 100
#'     x
#' }
#' ms <- matchSpectra(addProcessing(pest_ms2, norm_int),
#'                    addProcessing(minimb, norm_int), csp)
#'
#' ## validate matches using the shiny app. Note: the call is only executed
#' ## in interactive mode.
#' if (interactive()) {
#'     validateMatchedSpectra(ms)
#' }
validateMatchedSpectra <- function(object) {
    if (!requireNamespace("shiny", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'shiny'.",
             " Please install with 'BiocInstaller::install(\"shiny\")'")
    if (!requireNamespace("shinyjs", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'shinyjs'.",
             " Please install with 'BiocInstaller::install(\"shinyjs\")'")
    if (!requireNamespace("DT", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'DT'.",
             " Please install with 'BiocInstaller::install(\"DT\")'")
    if (!requireNamespace("SpectraVis", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'SpectraVis'.",
             " Please install with 'BiocInstaller::install(\"RforMassSpectrometry/SpectraVis\")'")

    stopifnot(inherits(object, "MatchedSpectra"))
    if (!length(object))
        stop("The 'MatchedSpectra' object is empty.")

    bpp <- bpparam()
    on.exit(register(bpp))
    register(SerialParam())
    ppm <- 20
    tolerance <- 0
    prm <- object@metadata$param
    if (length(prm)) {
        if (any(slotNames(prm) == "ppm"))
            ppm <- prm@ppm
        if (any(slotNames(prm) == "tolerance"))
            tolerance <- prm@tolerance
    }
    ui <- shiny::fluidPage(shinyjs::useShinyjs(),
                    shiny::titlePanel(shiny::div("validateMatchedSpectra")),
                    shiny::sidebarLayout(shiny::sidebarPanel(
                        shiny::selectInput("query", "Query Spectra:",
                                    choices = list(),
                                    multiple = FALSE,
                                    selectize = FALSE,
                                    size = 25),
                        shiny::actionButton("b_store", "Save & Close")
                    ),
                    shiny::mainPanel(
                        plotly::plotlyOutput("plot"),
                        shiny::checkboxInput("valid", "Current match OK?", value = TRUE,
                                      width = NULL),
                        DT::DTOutput("targets")
                        )
                    ))
    server <- function(input, output, session) {
        query_ids <- .createChoices(object)
        dt <- .create_dt(object)
        if (nrow(dt)) {
            dt <- data.frame(valid = TRUE, dt)
        } else dt <- data.frame(valid = logical(), dt)
        dtl <- split(dt, factor(object@matches$query_idx, seq_along(object)))
        rv <- shiny::reactiveValues(
            queries = query_ids,
            dtl = dtl
            )
        rv_query <- shiny::reactiveValues(idx = 1L)
        rv_target <- shiny::reactiveValues(idx = 1L)
        shiny::observe(shiny::updateSelectInput(inputId = "query",
                                  choices = rv$queries))
        shinyjs::disable("valid")
        ## Choose query spectrum
        shiny::observeEvent(input$query, {
            rv_query$idx <- as.integer(input$query)
            rv_target$idx <- 1L
            current_match <- object[rv_query$idx]
            output$targets <- DT::renderDT(
                DT::datatable(rv$dtl[[rv_query$idx]],
                          selection = list(mode = "single",
                                           selected = rv_target$idx),
                          options = list(dom = "t", ordering = FALSE),
                          escape = FALSE, rownames = FALSE),
                server = TRUE)
            if (nrow(rv$dtl[[rv_query$idx]]) >= rv_target$idx) {
                ## rv_target$idx <- input$targets_rows_selected
                shinyjs::enable("valid")
                current_valid <- rv$dtl[[rv_query$idx]]$valid[rv_target$idx]
                shiny::updateCheckboxInput(session, "valid", value = current_valid)
            } else
                shinyjs::disable("valid")
            output$plot <- plotly::renderPlotly(
                                       SpectraVis::plotlySpectraMirror(
                                                       query(current_match),
                                                       Spectra(),
                                                       xLabel = "query",
                                                       yLabel = "target",
                                                       ppm = ppm,
                                                       tolerance = tolerance))
        })
        ## Choose target spectrum
        shiny::observeEvent(input$targets_rows_selected, {
            rv_target$idx <- input$targets_rows_selected
            current_match <- object[rv_query$idx]
            if (nrow(current_match@matches)) {
                tidx <- current_match@matches$target_idx[rv_target$idx]
                current_valid <- rv$dtl[[rv_query$idx]]$valid[rv_target$idx]
                shiny::updateCheckboxInput(session, "valid", value = current_valid)
                output$plot <- plotly::renderPlotly(
                                           SpectraVis::plotlySpectraMirror(
                                                           query(current_match),
                                                           target(current_match)[tidx],
                                                           xLabel = "query",
                                                           yLabel = "target",
                                                           ppm = ppm,
                                                           tolerance = tolerance))
            }
        })
        shiny::observeEvent(input$valid, {
            if (length(rv_target$idx) &&
                nrow(rv$dtl[[rv_query$idx]]) >= rv_target$idx) {
                rv$dtl[[rv_query$idx]]$valid[rv_target$idx] <- input$valid
            }
        })
        shiny::observeEvent(input$b_store, {
            idx <- which(do.call(rbind, rv$dtl)$valid)
            shiny::stopApp(filterMatches(object, index = idx))
        })
    }
    shiny::runApp(shiny::shinyApp(ui, server))
}

#' isolate entries for feature selection
#' @noRd
.createChoices <- function(x) {
    sl <- seq_along(x)
    l_choices <- as.list(sl)
    names(l_choices) <- paste0(sl, " - MZ",
                               round(query(x)$precursorMz, 4), "@RT",
                               round(query(x)$rtime / 60, 2), " min")
    l_choices
}

.create_dt <- function(x){
    .sel_cols <- c("precursorMz", "target_precursorMz", "rtime",
                   "target_rtime", "target_name", "target_compound_name",
                   "score", "reverse_score", "presence_ratio")
    cols <- .sel_cols[.sel_cols %in% spectraVariables(x)]
    tbl <- as.data.frame(matchedData(x, cols))
    tbl$score <- round(tbl$score, 3)
    if (any(colnames(tbl) == "reverse_score"))
        tbl$reverse_score <- round(tbl$reverse_score, 3)
    if (any(colnames(tbl) == "presence_ratio"))
        tbl$presence_ratio <- round(tbl$presence_ratio, 3)
    if (any(colnames(tbl) == "rtime"))
        tbl$rtime <- round(tbl$rtime, 1)
    if (any(colnames(tbl) == "target_rtime"))
        tbl$target_rtime <- round(tbl$target_rtime, 1)
    tbl[!is.na(tbl$score), , drop = FALSE]
}