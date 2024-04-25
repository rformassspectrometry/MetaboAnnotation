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
#' Note that column `"query_index_"` and `"target_index_"` are temporarily
#' added to the query and target `Spectra` object to display them in the
#' interactive graphics for easier identification of the compared spectra.
#'
#' @param object A non-empty instance of class `MatchedSpectra`.
#'
#' @return A `MatchedSpectra` with validated results.
#'
#' @export
#'
#' @importFrom BiocParallel bpparam SerialParam register
#'
#' @importFrom methods slotNames
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
#'     res <- validateMatchedSpectra(ms)
#' }
validateMatchedSpectra <- function(object) {
    if (!requireNamespace("shiny", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'shiny'.",
             " Please install with 'BiocManager::install(\"shiny\")'")
    if (!requireNamespace("shinyjs", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'shinyjs'.",
             " Please install with 'BiocManager::install(\"shinyjs\")'")
    if (!requireNamespace("DT", quietly = TRUE))
        stop("The use of 'validateMatchedSpectra' requires package 'DT'.",
             " Please install with 'BiocManager::install(\"DT\")'")

    stopifnot(inherits(object, "MatchedSpectra"))
    if (!length(object))
        stop("The 'MatchedSpectra' object is empty.")

    ## Add query and target index
    object@query$query_index_ <- seq_along(object@query)
    object@target$index_ <- seq_along(object@target)

    bpp <- bpparam()
    on.exit(register(bpp))
    register(SerialParam())
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
                        shiny::checkboxInput("valid", "Current match OK?",
                                             value = TRUE, width = NULL),
                        DT::DTOutput("targets")
                        )
                    ))
    server <- function(input, output, session) {
        query_ids <- .createChoices(object)
        dt <- .create_dt(object)
        ppm <- 20
        tolerance <- 0
        if (length(prm <- object@metadata$param)) {
            if ("ppm" %in% slotNames(prm))
                ppm <- prm@ppm
            if ("tolerance" %in% slotNames(prm))
                tolerance <- prm@tolerance
        }
        if (nrow(dt)) {
            dt <- data.frame(valid = TRUE, dt)
        } else dt <- data.frame(valid = logical(), dt)
        dtl <- split(dt, factor(dt$query_index_, seq_along(object)))
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
                shiny::updateCheckboxInput(session, "valid",
                                           value = current_valid)
            } else
                shinyjs::disable("valid")
            output$plot <- plotly::renderPlotly(
                        .plotlySpectraMirror(query(current_match), Spectra(),
                                             xLabel = "query",
                                             xColor = "#E41A1C",
                                             yLabel = "target",
                                             yColor = "#377EB8"))
        })
        ## Choose target spectrum
        shiny::observeEvent(input$targets_rows_selected, {
            rv_target$idx <- input$targets_rows_selected
            current_match <- object[rv_query$idx]
            if (nrow(current_match@matches)) {
                tidx <- current_match@matches$target_idx[rv_target$idx]
                current_valid <- rv$dtl[[rv_query$idx]]$valid[rv_target$idx]
                shiny::updateCheckboxInput(session, "valid",
                                           value = current_valid)
                output$plot <- plotly::renderPlotly(
                        .plotlySpectraMirror(query(current_match),
                                             target(current_match)[tidx],
                                             ppm = ppm, tolerance = tolerance,
                                             xLabel = "query",
                                             xColor = "#E41A1C",
                                             yLabel = "target",
                                             yColor = "#377EB8"))
            }
        })
        shiny::observeEvent(input$valid, {
            if (length(rv_target$idx) &&
                nrow(rv$dtl[[rv_query$idx]]) >= rv_target$idx) {
                rv$dtl[[rv_query$idx]]$valid[rv_target$idx] <- input$valid
            }
        })
        shiny::observeEvent(input$b_store, {
            ## Collect all the selections from all data frames
            idx <- which(do.call(rbind, rv$dtl)$valid)
            shiny::stopApp(filterMatches(object, index = idx))
        })
    }
    shiny::runApp(shiny::shinyApp(ui, server))
}

#' function to create interactive plot
#'
#' This is `SpectraVis::plotlySpectraMirror` - replace with import or require
#' once `SpectraVis` is in Bioconductor.
#'
#' @importMethodsFrom Spectra peaksData
#'
#' @importFrom MsCoreUtils common
#'
#' @noRd
.plotlySpectraMirror <- function(x, y, xLabel = "", xColor = "#737373",
                                yLabel = "", yColor = "#737373", matchSize = 5,
                                ppm = 20, tolerance = 0) {
    stopifnot(inherits(x, "Spectra"))
    stopifnot(inherits(y, "Spectra"))
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("The use of '.plotlySpectraMirror' requires package 'plotly'. ",
             "Please install with 'BiocInstaller::install(\"plotly\")'")
    p <- plotly::plot_ly()
    if (length(x) > 1 || length(y) > 1)
        stop("'x' and 'y' have to be of length 1")
    if (length(x))
        x_peaks <- as.data.frame(peaksData(x)[[1L]])
    else x_peaks <- data.frame(mz = numeric(), intensity = numeric(),
                               match = character())
    if (length(y))
        y_peaks <- as.data.frame(peaksData(y)[[1L]])
    else y_peaks <- data.frame(mz = numeric(), intensity = numeric(),
                               match = character())

    x_range <- range(x_peaks$mz, y_peaks$mz, na.rm = TRUE) + c(-1, 1)
    y_max <- max(x_peaks$intensity, y_peaks$intensity, na.rm = TRUE)
    y_peaks$intensity <- -y_peaks$intensity

    ht <- "<b>%{text}</b><br>mz: %{x}<br>int: %{y}"
    if (nrow(x_peaks)) {
        x_peaks$zero <- 0.0
        x_peaks$match <- ""
        x_peaks$color <- xColor[1L]
        idx <- which(common(x_peaks$mz, y_peaks$mz, tolerance, ppm))
        if (length(idx))
            x_peaks$match[idx] <- "matched"
        p <- .plotly_peaks(p, x_peaks, name = xLabel, col = xColor[1L],
                           hovertemplate = ht, text = ~match)
    }
    if (nrow(y_peaks)) {
        y_peaks$zero <- 0.0
        y_peaks$match <- ""
        y_peaks$color <- yColor[1L]
        idx <- which(common(y_peaks$mz, x_peaks$mz, tolerance, ppm))
        if (length(idx))
            y_peaks$match[idx] <- "matched"
        p <- .plotly_peaks(p, y_peaks, name = yLabel, col = yColor[1L],
                           hovertemplate = ht, text = ~match)
    }
    pks <- rbind(x_peaks, y_peaks)
    pks <- pks[pks$match != "", , drop = FALSE]
    if (nrow(pks))
        p <- plotly::add_trace(p, data = pks, x = ~mz, y = ~intensity,
                       type = "scatter", mode = "markers",
                       hoverinfo = "none", name = "matched",
                       marker = list(size = matchSize[1L],
                                     color = ~color))
    plotly::layout(p, xaxis = list(title = "m/z", zeroline = FALSE),
           yaxis = list(title = "intensity", zeroline = TRUE),
           hovermode = "x", hoverdistance = 1)
}

#' That's also from `SpectraVis` - remove if `SpectraVis` is used instead.
#'
#' @noRd
.plotly_peaks <- function(p, data, col = "#737373", name = "",
                          hovertemplate = "<br>mz: %{x}<br>int: %{y}<br>",
                          ...) {
    plotly::add_segments(p, data = data, x = ~mz, y = ~zero, xend = ~mz,
                 yend = ~intensity, line = list(color = col),
                 name = name, hovertemplate = hovertemplate, ...)
}


## .plotlySpectraMirror <- function(x, y, col = c("#E41A1C", "#377EB8"),
##                                  main = "") {
##     if (!requireNamespace("plotly", quietly = TRUE))
##         stop("The use of '.plotlySpectraMirror' requires package 'plotly'. ",
##              "Please install with 'BiocInstaller::install(\"plotly\")'")
##     if (length(col) != 2)
##         col <- col[c(1, 1)]
##     if (length(x)) {
##         upper <- as.data.frame(peaksData(x)[[1L]])
##     } else upper <- data.frame(mz = numeric(), intensity = numeric())
##     if (length(y)) {
##         lower <- as.data.frame(peaksData(y)[[1L]])
##     } else lower <- data.frame(mz = numeric(), intensity = numeric())
##     if (nrow(upper))
##         upper$zero <- 0.0
##     else upper$zero <- numeric()
##     if (nrow(lower))
##         lower$zero <- 0.0
##     else lower$zero <- numeric()

##     mz_range <- range(upper$mz, lower$mz) + c(-1, 1)
##     maxy <- max(upper$intensity, lower$intensity, na.rm = TRUE)
##     int_range <- list(-maxy, maxy)

##     lower$intensity <- -lower$intensity
##     p <- plotly::plot_ly()
##     p <- plotly::add_segments(p, data = upper, x = ~mz, y = ~zero, xend = ~mz,
##                               yend = ~intensity, line = list(color = col[1L]),
##                               name = "query",
##                               hovertemplate = "<br>mz: %{x}<br>int: %{y}<br>")
##     p <- plotly::add_segments(p, data = lower, x = ~mz, y = ~zero, xend = ~mz,
##                               yend = ~intensity, line = list(color = col[2L]),
##                               name = "target",
##                               hovertemplate = "<br>mz: %{x}<br>int: %{y}<br>")
##     p <- plotly::layout(p, title = main,
##                         xaxis = list(title = "m/z", range = mz_range),
##                         yaxis = list(title = "intensity", zeroline = TRUE,
##                                      range = int_range),
##                         hovermode = "x", hoverdistance = 1)
##     p
## }

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
    .sel_cols <- c("query_index_", "target_index_", "precursorMz",
                   "target_precursorMz", "rtime", "target_rtime",
                   "target_name", "target_compound_name",
                   "score", "reverse_score", "presence_ratio")
    cols <- .sel_cols[.sel_cols %in% spectraVariables(x)]
    tbl <- as.data.frame(matchedData(x, cols))
    tbl$score <- round(tbl$score, 3)
    if (any(colnames(tbl) == "precursorMz"))
        tbl$precursorMz <- round(tbl$precursorMz, 3)
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
