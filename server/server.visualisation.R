
  # Creating and populating the x axis
  output$x.axis <- renderUI({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    selectInput("x.axis", "X axis",names(df))
  })
  
  # Creating and populating the y axis
  output$y.axis <- renderUI({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    selectInput("y.axis", "Y axis",names(df))
  })
  
  # Creating and populating the statistics axis
  output$summary.axis <- renderUI({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    selectInput("summary.axis", "Select a variable for statistics",names(df))
  })
  
  # Scatter plot
  output$scatterPlot <- renderPlotly({
    # features to be plotted
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    fig <- plot_ly(data=df,
                   x=df[,input$x.axis],
                   y=df[,input$y.axis], type = 'scatter', mode='markers')
    fig <- fig %>% layout(xaxis = list(title = input$x.axis),
                          yaxis = list(title = input$y.axis))
  })
  
  # Density plot
  output$densityPlot <- renderPlotly({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    fig <- plot_ly(x=df[,input$x.axis],
                   y=df[,input$y.axis]) %>% 
      add_histogram2dcontour(showscale=TRUE, ncontours=20, colorscale='Viridis', 
                             contours = list(coloring='heatmap')) %>%
      add_markers(x = df[,input$x.axis], 
                  y = df[,input$y.axis], marker=list(size=3))
    fig <- fig %>% layout(xaxis = list(title = input$x.axis),
                          yaxis = list(title = input$y.axis))
  })
  
  # Data summary text
  output$summaryTab <- renderPrint({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    summary(df[,input$summary.axis])
  })
  
  # Boxplot
  output$boxplot <- renderPlotly({
    df <- rv$data.sample$features$preprocessed$x
    if (is.null(df)) return(NULL)
    fig <- plot_ly(y = df[,input$summary.axis], type = "box", name = input$summary.axis)
  })
  
  # Signals slider
  output$sig.slider <- renderUI({
    df_sig_1 <<- rv$data.sample$profiles
    # Converting from "by" class 
    df_sig <-lapply(df_sig_1,function(x){ e<-rbind.fill(x)})
    if (is.null(df)) return(NULL)
    View(df_sig)
    sliderInput("sig.slider", "Observations",min = 0, max = length(names(df_sig)), value = 0, step = 1, width = '100%')
  })
  
  # Signals plotting
  observeEvent(input$sig.slider,
     output$signals <- renderPlotly({
       df_sig_1 <- rv$data.sample$profiles
       # Converting from "by" class 
       df_sig <-lapply(df_sig_1,function(x){ e<-rbind.fill(x)})
       
       if (is.null(df_sig)) return(NULL)
       data_plot(df_sig,as.character(input$sig.slider))    
               })        
               )
  
  
  