  
  
  tabPanel("Visualisation",
    sidebarLayout(
      # Selectors
      sidebarPanel(
        uiOutput("x.axis"),
        uiOutput("y.axis"),
        uiOutput("summary.axis")
      ),
      mainPanel(
        # Plots
        tabsetPanel(
          tabPanel("Scatter Plot",
                   plotlyOutput("scatterPlot")
                   ),
          tabPanel("Density Plot",
                   plotlyOutput("densityPlot")
                   ),
          tabPanel("Variable statistics",
                   splitLayout(
                     plotlyOutput("boxplot"),
                     verbatimTextOutput("summaryTab")
                   )),
          tabPanel("Signals",
                   plotlyOutput("signals"),
                   uiOutput("sig.slider"))
        )
      )
    )
  )