  
  tabPanel(title = "Importation",
    # Messages collapse panel
    sidebarLayout(
      sidebarPanel(width = 3, 
        bsCollapsePanel("Messages", verbatimTextOutput("messages"), style = "primary")
      ),
      
      mainPanel(
        position = 'left',
        # Loading the csv features file
        tags$h4("REQUIRED INPUT DATA FILES"),
        hr(),
        fluidRow(
          column(width = 4,
                 p(HTML("<b>FEATURES/RDS FILE</b>"),span(shiny::icon("info-circle"), id = "csvHelp"),
                   tippy::tippy_this(elementId = "csvHelp",tooltip = "<span style='font-size:15px;'>Data must be in a .csv file. Observations in rows. Features in columns. Missing value must be empty. After a first use of your files an Rclustool RDS file is saved. This file contains all the data you used and it's faster to load",placement = "right")
                 )
          )),
        fluidRow(
          column(width = 2,
                 shinyFilesButton("featuresButton", "FEATURES/RDS FILE", "Please select a file", 
                                  multiple = FALSE, viewtype = "detail", icon = icon("file-csv"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("featuresText")
          ),
          column(width = 1, offset = 1,
                 selectInput("sepFeat","Sep",choices =c(",",";","\\s","\\t"), selected = ",", width = '65%')
          ),
          column(width = 1,
                 selectInput("decFeat","Dec",choices =c(".",","), selected = ".", width = '65%')
          ),
          column(width = 1,
                 selectInput("naFeat","Missing",choices =c(" ","NA","9999"), selected = " ", width = '65%')
          )
        ),
        
        
        # Loading optional data files
        tags$h4("OPTIONAL INPUT DATA FILES"),
        hr(),
        fluidRow(
          column(width = 4,
                 p(HTML("<b>SIGNALS</b>"),span(shiny::icon("info-circle"), id = "signalsHelp"),
                   tippy::tippy_this(elementId = "signalsHelp",tooltip = "<span style='font-size:15px;'> Data must be in a .csv file. Signals in columns. A same ID for all signal values. Missing value must be empty.",placement = "right")
                 )
          )),
        # Signals
        fluidRow(
          column(width = 2,
                 shinyFilesButton("signalsButton", "SIGNALS", "Please select a file", 
                                  multiple = FALSE, viewtype = "detail", icon = icon("file-csv"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("signalsText")
          ),
          column(width = 1, offset = 1,
                 selectInput("sepSig","Sep",choices =c(",",";","\\s","\\t"), selected = ",", width = '65%')
          ),
          column(width = 1,
                 selectInput("decSig","Dec",choices =c(".",","), selected = ".", width = '65%')
          ),
          column(width = 1,
                 selectInput("naSig","Missing",choices =c(" ","NA","9999"), selected = " ", width = '65%')
          )
        ),
        hr(),
        
        # Images
        fluidRow(
          column(width = 4,
                 p(HTML("<b>IMAGES</b>"),span(shiny::icon("info-circle"), id = "imagesHelp"),
                   tippy::tippy_this(elementId = "imagesHelp",tooltip = "<span style='font-size:15px;'> JPEG or PNG images. Observation's ID for file name.",placement = "right")
                 )
          )),
        fluidRow(
          column(width = 2,
                 actionButton(inputId = "imageButton", label = "IMAGES",
                              icon = icon("folder"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("imageDirtext"))
        ),
        hr(),
        
        # Metadata
        fluidRow(
          column(width = 4,
                 p(HTML("<b>METADATA</b>"),span(shiny::icon("info-circle"), id = "metadataHelp"),
                   tippy::tippy_this(elementId = "metadataHelp",tooltip ="<span style='font-size:15px;'> Data must be in a .txt file. \"Metadata name: value\"",placement = "right")
                 )
          )),
        fluidRow(
          column(width = 2,
                 shinyFilesButton("metadataButton", "METADATA", "Please select a file", 
                                  multiple = FALSE, viewtype = "detail", icon = icon("file-alt"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("metadataText"))
        ),
        hr(),
        
        # Working directory
        tags$h4("WORKING DIRECTORY"),
        hr(),
        fluidRow(
          column(width = 2,
                 actionButton(inputId = "setwdButton", label = "DIRECTORY",
                              icon = icon("folder"))
          ),
          column(width = 4, offset = 1,
                 verbatimTextOutput("wdText"))
        ),        
        hr(),
        # Import and Reset Buttons
        fluidRow(
          column(width = 2, 
                 actionButton(inputId = "importButton", label = "IMPORT",width = "95px", icon = icon("database"))
          ),
          column(width = 2, offset = 2, 
                 actionButton(inputId = "resetButton", label = "RESET",width = "95px", icon = icon("undo"))
          )
        ),
        hr()
      ) # MainPanel bracket
    )
              
            
  
  ) #tabPanel bracket
  
  
  
