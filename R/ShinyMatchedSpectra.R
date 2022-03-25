#' @title Verify Results from matchSpectra
#'
#' @description
#'
#' The `shinyMatchedSpectra()` function opens a simple shiny application
#' that allows to browse results stored in a `MatchedSpectra` object. Results
#' can be verified and set to TRUE or FALSE. Upon pushing the "Save & Close"
#' button the the app is closed and a filtered `MatchedSpectra` is returned,
#' containing only results set to TRUE.
#'
#' @param object A non-empty instance of class `MatchedSpectra`.
#'
#' @return A `MatchedSpectra` with verified results
#'
#' @export
#'
#' @import shiny
#' @import DT
#' @import BiocParallel
#' @import Spectra
#' @import S4Vectors
#'
#' @author Carolin Huber, Michael Witting
shinyMatchedSpectra <- function(object) {

  stopifnot(inherits(object, "MatchedSpectra"))
  if (!length(object)) {
    stop("The 'MatchedSpectra' object is empty.")
  }

  # set serial processing (important for fast shiny app)
  register(SerialParam())


  ui <- fluidPage(
    titlePanel(
      div("ShinyMatchedSpectra")
    ),
    sidebarLayout(
      sidebarPanel(
        selectInput("selection", "Query Spectra:",
                    choices = list(),
                    multiple = FALSE,
                    selectize = FALSE,
                    size = 25),
        actionButton("b_store", "Save & Close")
      ),
      mainPanel(
        plotlyOutput("plot"),
        DT::dataTableOutput("dynamic"),
        radioButtons("veri",
                     label = "Hit correct?",
                     choices = list("Yes" = T,
                                    "No" = F),
                     selected = T)
      )
    )
  )


  server <- function(input, output, session) {

    # reduce only to results with matches
    mtch_sub <- object[whichQuery(object)]
    mtch_sub <- pruneTarget(mtch_sub)

    q <- reactiveValues(mtch_sub = mtch_sub,
                        l_choices = .createChoices(mtch_sub),
                        boolean_values = .createBoolean(mtch_sub))

    observe(updateSelectInput(inputId='selection', choices = q$l_choices))

    # define storage places for reactive values that change during the shiny application run
    reacVal_mtch <- reactiveValues(mtch_sub = mtch_sub,
                              l_choices = .createChoices(mtch_sub),
                              boolean_values = .createBoolean(mtch_sub))

    reacVal_query <- reactiveValues(id = 1)
    reacVal_target <- reactiveValues(id = 1)

    specs <- reactiveValues(match = MatchedSpectra(),
                            query = Spectra(),
                            target = Spectra())

    observe(updateSelectInput(inputId = "selection",
                              choices = reacVal_mtch$l_choices))

    # Define what happens if a selection in list clicked
    observeEvent(input$selection, {

      # change index to selection
      reacVal_query$id <- as.numeric(input$selection)

      # filter MatchedSpectra based on index from selection
      x <- reacVal_mtch$mtch_sub[reacVal_query$id]
      mtch_tbl <- .createTable(x)

      # create dynamic table
      output$dynamic <- DT::renderDT(DT::datatable(as.data.frame(cbind.DataFrame(mtch_tbl,
                                                                                 Hit = reacVal_mtch$boolean_values[[reacVal_query$id]])),
                                                   #selection = "single"),
                                                   selection = list(mode = "single", selected = c(1))),
                                     server = TRUE)

      reacVal_target$id <- input$dynamic_rows_selected
      output$row <- renderPrint(input$dynamic_rows_selected)

      # update reactiveValues containing spectra
      specs$match <- x
      specs$query <- query(x)
      specs$target <- Spectra()

      # render new plotly plot
      observe(output$plot <-  renderPlotly(.plotlyMirrorPlot(specs$query, specs$target)))

      # change the button to the previous selected verification value
      updateRadioButtons(inputId="veri",
                         choices = list("Yes" = T,
                                        "No" = F),
                         selected = reacVal_mtch$boolean_values[reacVal_query$id][reacVal_target$id])
    })

    # Define what happens if a row in DT is selected
    observeEvent(input$dynamic_rows_selected, {

      # change index to selection
      reacVal_target$id <- input$dynamic_rows_selected

      # filter MatchedSpectra based on index from selection
      reacVal_query$id <- as.numeric(input$selection)
      x <- reacVal_mtch$mtch_sub[reacVal_query$id]

      # update reactiveValues containing spectra
      specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]

      # render new plotly plot
      observe(output$plot <- renderPlotly(.plotlyMirrorPlot(specs$query,
                                                            specs$target)))
    })

    # update verification
    observeEvent(input$veri, {
      reacVal_mtch$boolean_values[[reacVal_query$id]][reacVal_target$id] <- as.logical(input$veri)
    })

    # store results and close app
    observeEvent(input$b_store, {
      mtch_sub_filtered <- .filterVerified(reacVal_mtch$mtch_sub,
                                           reacVal_mtch$boolean_values)
      stopApp(mtch_sub_filtered)
    })

  }

  # runApp is required to make return of values working
  runApp(shinyApp(ui, server))

}

#' function to create interactive plot
#' @import plotly
#'
#' @noRd
.plotlyMirrorPlot <- function(query_spectrum, target_spectrum){

  # create data frames for plotting
  top <- data.frame(mz = unlist(mz(query_spectrum)),
                    int = unlist(intensity(query_spectrum)))

  top$int <- top$int / max(top$int) * 100

  if(length(target_spectrum) > 0) {

    bottom <- data.frame(mz = unlist(mz(target_spectrum)),
                         int = unlist(intensity(target_spectrum)))

    bottom$int <- bottom$int / max(bottom$int) * 100

    min_mz <- min(c(top$mz, bottom$mz))
    max_mz <- max(c(top$mz, bottom$mz))

  } else {

    min_mz <- min(top$mz)
    max_mz <- max(top$mz)

  }

  # create layout
  layout <- list(
    title = "",
    xaxis = list(title = "m/z",
                 zeroline = TRUE,
                 range = c(min_mz, max_mz),
                 nticks = 8,
                 autorange = TRUE),
    yaxis = list(title = "Signal Intensity [%]",
                 zeroline = TRUE,
                 tickmode ='array',
                 tickvals = c(-100, -50, 0, 50, 100),
                 ticktext = c('100','50', '0', '50', '100'))
  )

  # create plot
  p <- plot_ly(
    top,
    x =  ~ mz,
    y =  ~ int,
    showlegend = F,
    type = 'bar',
    marker = list(size = 3, color = 'red'),
    hoverinfo = 'none'
  )

  p <-
    add_markers(
      p,
      type = "scatter",
      x = top$mz,
      y = top$int,
      hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'),
      hoverlabel = list(namelength = 0)
    )

  if(length(target_spectrum) > 0) {

    p <-
      add_trace(
        p,
        type = "bar",
        x = bottom$mz,
        y = -bottom$int,
        marker = list(color = 'blue'),
        hoverinfo = 'none'
      )

    p <-
      add_markers(
        p,
        x = bottom$mz,
        y = -bottom$int,
        type = 'scatter',
        marker = list(color = 'blue'),
        hovertemplate = paste('<br>mz:', '%{x}', '<br>int: %{y}<br>'),
        hoverlabel = list(namelength = 0)
      )
  }

  p <-
    layout(
      p,
      title = layout$title,
      xaxis = layout$xaxis,
      yaxis = layout$yaxis
    )

  p <-
    add_annotations(
      p,
      type = 'text',
      x = c(15, 15),
      y = c(100, -100),
      text = c("query", "library"),
      textfont = list(color = c('red', 'blue')),
      showarrow = F
    )

  p <- p %>% layout(hovermode = "x", hoverdistance = 1)

  # return plot
  p
}

#' isolate entries for feature selection
#' @noRd
.createChoices <- function(mtch_sub){
  #generate choices name for list on side
  l_choices <- as.list(seq(1, length(mtch_sub)))
  names(l_choices) <- paste0(seq(1, length(mtch_sub)),
                             " - MZ",
                             round(query(mtch_sub)$precursorMz,4),
                             "@RT",
                             round(query(mtch_sub)$rtime/60, 2), " min")
  l_choices
}

#' create a list with vector of boolean values for verification
#' @noRd
.createBoolean <- function(mtch_sub){
  boolean_values <- list()
  for(i in 1:length(mtch_sub)) {
    x <- mtch_sub[i]
    target_length <- length(x@target[x@matches$target_idx])
    boolean_values[[i]] <- rep(TRUE, length(x@target[x@matches$target_idx]))
  }
  boolean_values
}

#' create table for display
#' @noRd
.createTable <- function(x){
    .sel_cols <- c("precursorMz", "target_precursorMz", "rtime",
                   "target_rtime", "target_name", "target_compound_name",
                   "score", "reverse_score", "presence_ratio")
    cols <- .sel_cols[.sel_cols %in% spectraVariables(x)]
    tbl <- matchedData(x, cols)

    tbl$score <- round(tbl$score, 3)
    tbl$reverse_score <- round(tbl$reverse_score, 3)
    tbl$presence_ratio <- round(tbl$presence_ratio, 3)

    tbl
}

#' filter matches based on manual selection
#' @noRd
.filterVerified <- function(mtch_sub, boolean_values) {

  filter_idx <- mtch_sub@matches[which(unlist(boolean_values)),]
  filterMatches(mtch_sub,
                index = as.integer(row.names(filter_idx)))
}
