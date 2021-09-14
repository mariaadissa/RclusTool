
  # Import tab environment
  prepro.env <<- new.env()
  RclusTool.env$gui$tabs.env$prepro <- c(prepro.env)
  
  # Testing preprocessText
  output$preprocessText <- renderText({
    paste("No file selected")
  })
  
  # Dynamically inserting visualization tab
  observeEvent(input$visualizeButton, {
    if(input$visualizeButton[1]==1){
      if(!is.null(rv$data.sample)){
        insertTab(inputId = "Rclustools",
                  source(file.path("ui", "ui.visualisation.R"), local = TRUE)$value,
                  target = "Preprocessing",
                  position = "after"
        )
      }
    }else{
      isolate(input$visualizeButton)
    }
  })