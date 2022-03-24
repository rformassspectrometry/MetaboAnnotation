#' @title Verify Results from matchSpectra
#'
#' @description
#'
#' The `shinyMatchedSpectra()` function opens a simple shiny application
#' that allows to browse
#'
#'
#' @param object A non-empty instance of class `MatchedSpectra`.
#'
#' @return A `MatchedSpectra` with verified results
#'
#' @export
#'
#' @import shiny
#' @import DT
#'
#' @author Carolin Huber, Michael Witting
shinyMatchedSpectra <- function(object) {

  stopifnot(inherits(object, "MatchedSpectra"))
  if (!length(object)) {
    stop("The 'MatchedSpectra' object is empty.")
  }

  # set serial processing (important for fast shiny app)
  register(SerialParam())


  # Define UI for application
  ui <- fluidPage(
    titlePanel(
      # title with logos
      div("ShinyMetaboAnnotation")
    ),
    sidebarLayout(
      # define sidebar with list of matches and true match false match buttons
      sidebarPanel(
        # define File upload
        # define Feature selection
        selectInput('selection', 'Features:', choices=list(), multiple=FALSE, selectize=FALSE, size=25),
        actionButton("b_store", "Close")
      ),
      # define main window with text, plotly plot
      mainPanel(
        plotlyOutput("plot"),
        DT::dataTableOutput("dynamic"),
        radioButtons("veri",
                     label = "Hit correct?",
                     choices = list("Yes" = T, "No" = F),
                     selected = T)
      )
    )
  )

  #######################################################################
  # Define server logic
  server <- function(input, output, session) {

    # reduce only to results with matches
    mtch_sub <- mtch[whichQuery(object)]
    mtch_sub <- pruneTarget(mtch_sub)

    q <- reactiveValues(mtch_sub = mtch_sub,
                        l_choices = .createChoices(mtch_sub),
                        boolean_values = .createBoolean(mtch_sub))


    observe(updateSelectInput(inputId='selection', choices = q$l_choices))


    # define storage places for reactive values that change during the shiny application run
    v <- reactiveValues(i = 1)
    id <- reactiveValues(target_idx = 1L)
    specs <- reactiveValues(match = MatchedSpectra(),
                            query = Spectra(),
                            target = Spectra())

    # Define what happens if a selection in list clicked
    observeEvent(input$selection, {

      # change index to selection
      v$i <- as.numeric(input$selection)

      # filter MatchedSpectra based on index from selection
      x <- q$mtch_sub[v$i]
      mtch_tbl <- .createTable(x)

      # create dynamic table
      output$dynamic <- DT::renderDT(DT::datatable(as.data.frame(cbind.DataFrame(mtch_tbl,
                                                                   Hit = q$boolean_values[[v$i]])),
                                                   selection = "single"),
                                     server = TRUE)

      id$target_idx <- input$dynamic_rows_selected
      output$row <- renderPrint(input$dynamic_rows_selected)

      # update reactiveValues containing spectra
      specs$match <- x
      specs$query <- query(x)
      specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]

      # render new plotly plot
      observe(output$plot <-  renderPlotly(.plotlyMirrorPlot(specs$query, specs$target)))

      # change the button to the previous selected verification value
      updateRadioButtons(inputId="veri",
                         choices = list("Yes" = T, "No" = F),
                         selected = q$boolean_values[v$i][id$target_idx])
    })

    # Define what happens if a row in DT is selected
    observeEvent(input$dynamic_rows_selected, {

      # change index to selection
      id$target_idx <- input$dynamic_rows_selected

      # filter MatchedSpectra based on index from selection
      v$i <- as.numeric(input$selection)
      x <- q$mtch_sub[v$i]

      # update reactiveValues containing spectra
      specs$target <- x@target[x@matches$target_idx][input$dynamic_rows_selected]

      # render new plotly plot
      observe(output$plot <-  renderPlotly(.plotlyMirrorPlot(specs$query, specs$target)))

      # change the button to the previous selected verification value
      updateRadioButtons(inputId="veri",
                         choices = list("Yes" = T, "No" = F),
                         selected = q$boolean_values[v$i][input$dynamic_rows_selected])

    })

    # Define what to do if click on radio button
    observeEvent(input$veri, {

      # store value in logical vector
      q$boolean_values[[v$i]][id$target_idx] <- as.logical(input$veri)

    })

    # Store Result from commandline
    observeEvent(input$b_store, {

      mtch_sub_filtered <- .filterVerified(q$mtch_sub, q$boolean_values)

      stopApp(mtch_sub_filtered)

    })




  }


  # runApp is required to make return of values working
  runApp(shinyApp(ui, server))

}

##' @import plotly
.plotlyMirrorPlot <- function(query_spectrum, target_spectrum){

  # create data frames for plotting
  top <- data.frame(mz = unlist(mz(query_spectrum)),
                    int = unlist(intensity(query_spectrum)))
  # create layout
  layout <- list(
    title = "",
    xaxis = list(title = "m/z",
                 zeroline=TRUE,
                 range=c(0, max(top$mz)),
                 nticks=8,
                 autorange = TRUE),
    yaxis = list(title = "Signal Intensity [%]",
                 zeroline=TRUE,
                 tickmode='array',
                 tickvals=c(-100, -50, 0, 50, 100),
                 ticktext=c('100','50', '0', '50', '100'))
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

    bottom <- data.frame(mz = unlist(mz(target_spectrum)),
                         int = unlist(intensity(target_spectrum)))

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

  tbl <- matchedData(x)
  tbl <- tbl[,c("precursorMz",
              "target_precursorMz",
              "rtime",
              "target_rtime",
              "target_name",
              "score",
              "reverse_score",
              "presence_ratio")]

  tbl$score <- round(tbl$score, 3)
  tbl$reverse_score <- round(tbl$reverse_score, 3)
  tbl$presence_ratio <- round(tbl$presence_ratio, 3)

  names(tbl)[names(tbl) == "precursorMz"] <- "Query m/z"
  names(tbl)[names(tbl) == "target_precursorMz"] <- "Target m/z"
  names(tbl)[names(tbl) == "rtime"] <- "Query RT"
  names(tbl)[names(tbl) == "target_rtime"] <- "Target RT"
  names(tbl)[names(tbl) == "target_name"] <- "Name"
  names(tbl)[names(tbl) == "score"] <- "Forward Score"
  names(tbl)[names(tbl) == "reverse_score"] <- "Reverse Score"
  names(tbl)[names(tbl) == "presence_ratio"] <- "Presence Ratio"

  tbl
}


.filterVerified <- function(mtch_sub, boolean_values) {

  filter_idx <- mtch_sub@matches[which(unlist(boolean_values)),]
  filterMatches(mtch_sub,
                index = as.integer(row.names(filter_idx)))
}
