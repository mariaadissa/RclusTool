  
  tabPanel(title = "Preprocessing",
    sidebarLayout(
      sidebarPanel(width = 3,
        bsCollapsePanel("Messages", verbatimTextOutput("msgs"), style = "primary")
      ),
      mainPanel(
        position = 'left',
        tags$h4("IMPORT PREPROCESSING / VISUALIZE DATA SAMPLE (OPTIONAL)"),
        hr(),
        # Import preprocessing
        fluidRow(
          column(width = 2,
                 shinyFilesButton("impPrepButton", "Import Preprocessing", "Please select a file", 
                                  multiple = FALSE, viewtype = "detail", icon = icon("file-csv"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("preprocessText")
          )
        ),
        # Visualization
        fluidRow(
          column(width = 2,
            actionButton("visualizeButton","Visualize Data Sample",icon = icon("bar-chart-o"))
          )
        ),
        # Variable selection
        tags$h4("VARIABLE SELECTION"),
        hr(),
        fluidRow(),
        fluidRow(
          column(width = 3,
            # these should be replaced by uiOutput elements and create them in the server using renderUI to populate the choices with the features
            selectInput("datasetVar","Dataset variables:", choices = list("V1" = 1, "V2" = 2, "V3" = 3), selected = 1)),
          column(width = 3, offset = 1,
            selectInput("classifVar","Variables used for classification",choices = list("V1" = 1, "V2" = 2, "V3" = 3), selected = 1))
        ),
        fluidRow(
          column(width = 1,offset=3,
                 actionButton("toButton","=>"))
        ),
        fluidRow(
          column(width = 1,offset=3,
                 actionButton("fromButton","<="))
        ),
        hr(),
        wellPanel(
                 actionButton("filterButton","Filter"),
                 actionButton("transformButton","Transform"),
                 actionButton("createButton","Create"),
                 actionButton("configButton","Configure default plot"), 
                 style = "padding: 5px;"
        ),
        hr(),
        fluidRow(
          column(width = 1,
                 actionButton("expPrepButton","Export Preprocessing",icon = icon("file-csv"))),
          column(width = 1,offset = 2,
                 actionButton("validatePrepButton","Validate Preprocessing",icon = icon("fas fa-check"))),
          column(width = 1,offset = 3,
                 actionButton("resetPrepButton","Reset",icon = icon("undo")))
        )
        
      ) #mainPanel closing bracket
        
    ) #sidebarLayour closing bracket
  )