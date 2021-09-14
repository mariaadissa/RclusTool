
library(shiny)

# UI libraries
library(shinythemes)
library(shinyBS)
library(tippy)
library(bslib)

# Server libraries
library(easycsv)
library(shinyFiles)
library(sjmisc)
library(plotly)
library(plyr)

#Unused
library(pryr)

#functions
source(file.path("function", "parameters.R"))
source(file.path("function", "sampleBuilding.R"))
source(file.path("function", "sampleClustering.R"))
source(file.path("function", "data_plot.R"))

#global
source(file.path("global.R"))

  ui <- list( ui = tagList(
    
    tags$link( rel = "stylesheet", type = "text/css", href = "styles.css")),
    div(
      id = "mainPage",
      navbarPage(
        theme = shinytheme("cerulean"),
        title='RclusTools',
        id="Rclustools",
        source(file.path("ui", "ui.importation.R"), local = TRUE)$value,
        source(file.path("ui", "ui.preprocessing.R"), local = TRUE)$value
      )
    ))
  
  server <- function(input, output, session){
    
    # Setting up RclusTool environment
    if(!length(RclusTool.env)){
      initParameters(RclusTool.env)
    }
    
    source(file.path("server", "server.importation.R"), local = TRUE)$value
    source(file.path("server", "server.preprocessing.R"), local = TRUE)$value
    source(file.path("server", "server.visualisation.R"), local = TRUE)$value
    
  }
  
  # Build the application
  shinyApp(ui = ui, server = server)

