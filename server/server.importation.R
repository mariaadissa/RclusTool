  
  # Reactive values list
  rv <<- reactiveValues(rdsFile="",featuresFile="",Dir="",DirDefault="",
                        signalFile="",metadataFile="",imageDir="",data.sample=NULL)
  
  # Import tab environment
  import.env <<- new.env()
  RclusTool.env$gui$tabs.env <- c(import.env)
  
  # getting user's root directory
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home())

  #### FEATURES IMPORT #### 
  observe({
    shinyFileChoose(input, "featuresButton", roots = volumes, session = session, filetypes=c('RDS', 'csv'))
    #Get path
    if(!is.null(input$featuresButton)){
      import.env$featuresFile <- parseFilePaths(volumes, input$featuresButton)
      # extracting the directory name from chosen file
      if(!is_empty(import.env$featuresFile)){
        import.env$rdsFile <- ""
        rv$rdsFile <- ""  #reactiveVal
        import.env$Dir <- dirname(as.data.frame(import.env$featuresFile)$datapath)
        rv$featuresFile <- import.env$featuresFile$datapath #reactiveVal
        rv$Dir <- dirname(as.data.frame(import.env$featuresFile)$datapath) #reactiveVal
        import.env$DirDefault <- import.env$Dir
        rv$DirDefault <- import.env$Dir #reactiveVal
        # if an rds file is loaded
        if(grepl('.RDS',as.character(import.env$featuresFile$datapath))){
          import.env$rdsFile <- as.data.frame(import.env$featuresFile)$datapath
          rv$rdsFile <- as.data.frame(import.env$featuresFile)$datapath #reactiveVal
          import.env$featuresFile$datapath <- ""
          rv$featuresFile <- "" #reactiveVal
        }
      }
      # printing in the text field
      output$featuresText <- renderText({
        if(is_empty(import.env$featuresFile)){
          if(is_empty(import.env$rdsFile)){
            paste("No file selected")
          }
        }else if(import.env$rdsFile!=""){
          paste(import.env$rdsFile)
        }else {
          as.character(import.env$featuresFile$datapath)
        }
      })
    }
  })
  
  
  
  #### Signals import ####
  observe({
    shinyFileChoose(input,"signalsButton", roots = volumes, session = session, filetypes='csv')
    #Get path 
    if(!is.null(input$signalsButton)){
      import.env$signalFile <- parseFilePaths(volumes, input$signalsButton)
      import.env$signalFile <- as.data.frame(import.env$signalFile)
      rv$signalFile <- import.env$signalFile$datapath #reactiveVal
      # Verification of the existence of a datapath
      if(identical(import.env$signalFile$datapath, character(0))){
        import.env$signalFile[1,"datapath"] <- ""
        rv$signalFile <- "" #reactiveVal
      }
    }
    # printing in the text field 
    output$signalsText <- renderText({
      if(import.env$signalFile$datapath == ""){
        paste("No file selected")
      }else{
        as.character(import.env$signalFile$datapath)
      }
    })
  })
  
  
  #### Metadata import ####
  observe({
    shinyFileChoose(input,"metadataButton", roots = volumes, session = session, filetypes='txt')
    #Get path
    if(!is.null(input$metadataButton)){
      import.env$metadataFile <- parseFilePaths(volumes, input$metadataButton)
      import.env$metadataFile <- as.data.frame(import.env$metadataFile)
      rv$metadataFile <- import.env$metadataFile$datapath #reactiveVal
      # Verification of the existence of a datapath
      if(identical(import.env$metadataFile$datapath, character(0))){
        import.env$metadataFile[1,"datapath"] <- ""
        rv$metadataFile <- "" #reactiveVal
      }
    }
    # printing in the text field 
    output$metadataText <- renderText({
      if(import.env$metadataFile$datapath == ""){
        paste("No file selected")
      }else{
        as.character(import.env$metadataFile$datapath)
      }
    })
  })
  
  
  #### Images import ####
  import.env$imageDir <- ""
  observeEvent(input$imageButton, {
    # chosen image directory
    imdir <- choose_dir()
    # printing in the text field
    output$imageDirtext <- renderText({
      if(exists("imdir") && !is.na(imdir)){
        paste(imdir)
        import.env$imageDir <- imdir
        rv$imageDir <- import.env$imageDir #reactiveVal
      } else{paste("No images directory selected!")}
    })
  })

  
  # Working Directory
  import.env$Dir <- NULL
  observeEvent(input$setwdButton, {
    # chosen image directory
    wdir <- choose_dir()
    import.env$Dir <- wdir
    rv$Dir <- import.env$Dir  #reactiveVal
    # printing in the text field
    output$wdText <- renderText({
      if(is_empty(wdir)){wdir <- NA}
      if(exists("wdir") && !is.na(wdir)){paste(wdir)} 
      else{paste("No working directory selected!")}
    })
  })
  
  # Messages section
  textMsgs <- reactiveVal("")
  output$messages <- renderText(
    textMsgs()
  )
  
  ### Import Button ###
  observeEvent(input$importButton,{
    # No features or rds file is loaded
    if(is_empty(import.env$featuresFile) && is_empty(import.env$rdsFile)){
      # Error messages pop up window
      showModal(modalDialog(
        title = "Error message",
        paste0("No data to import. \nplease select a features file!"),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if(sum(is.na(import.env$featuresFile))==0 || !is.na(import.env$rdsFile)){
      
      importSample(file.features = import.env$featuresFile$datapath, file.meta = import.env$metadataFile$datapath,
                   file.profiles = import.env$signalFile$datapath, dir.images = import.env$imageDir, dir.save = import.env$DirDefault,
                   file.RDS = import.env$rdsFile, sepFeat = input$sepFeat, decFeat = input$decFeat,
                   naFeat=input$naFeat, sepSig = input$sepSig, decSig= input$decSig,
                   naSig=input$naSig, RclusTool.env=RclusTool.env) -> RclusTool.env$data.sample
      #reactiveVal
      #rv$data.sample <- RclusTool.env$data.sample
      importSample(file.features = rv$featuresFile, file.meta = rv$metadataFile,
                   file.profiles = rv$signalFile, dir.images = rv$imageDir, dir.save = rv$DirDefault,
                   file.RDS = rv$rdsFile, sepFeat = input$sepFeat, decFeat = input$decFeat,
                   naFeat=input$naFeat, sepSig = input$sepSig, decSig= input$decSig,
                   naSig=input$naSig, RclusTool.env=RclusTool.env) -> rv$data.sample
      
      # Displaying information in messages section
      if (!is.null(RclusTool.env$data.sample)){
        textMsgs(paste(textMsgs(),paste("----- Features importation -----\n",
                                        "Filename:  ", import.env$featuresFile$name, "\n",
                                        "Number of observations:  ", RclusTool.env$data.sample$size, "\n",
                                        "Number of features:  ", length(RclusTool.env$data.sample$features$initial$x), "\n\n", sep = "")))
        if (nchar(import.env$metadataFile$datapath)) {
          textMsgs(paste(textMsgs(), paste("----- MetaData importation -----\n",
                                           "Filename:  ", basename(RclusTool.env$data.sample$files$meta), "\n", sep = "")))
        }
        if (nchar(import.env$signalFile$datapath)) {
          textMsgs(paste(textMsgs(), paste("----- Signals importation -----\n",
                                           "Filename:  ", basename(RclusTool.env$data.sample$files$profiles), "\n",
                                           "Number of observations:  ", length(RclusTool.env$data.sample$profiles), "\n\n", sep = "")))
        }
        if (nchar(import.env$imageDir)) {
          textMsgs(paste(textMsgs(), paste("----- Images importation -----\n",
                                           "Folder:  ", basename(RclusTool.env$data.sample$files$images), "\n",
                                           "Number of observations:  ", length(!is.na(RclusTool.env$data.sample$images)), "\n\n", sep = "")))
        }
        if (nchar(import.env$DirDefault)) {
          textMsgs(paste(textMsgs(), paste("----- Working directory -----\n",
                                           "Folder:  ", basename(RclusTool.env$data.sample$files$dir), "\n\n")))
        }
        if (nchar(import.env$rdsFile)) {
          textMsgs(paste(textMsgs(), paste("----- RDS file -----\n",
                                           "Creation of a RDS object\n\n", sep = "")))
        }
      }
    }
  })# import observeEvent closing bracket
  
  
  ### Reset Button ###
  observeEvent(input$resetButton,{
    
    rv$data.sample <- NULL #reactiveVal
    
    # resetting features file datapath and text field
    import.env$featuresFile$datapath <- ""
    output$featuresText <- renderText({
          paste("No file selected")
    })
    
    # resetting signals file datapath and text field
    import.env$signalFile$datapath <- ""
    output$signalsText <- renderText({
      paste("No file selected")
    })
    
    # resetting metadata file datapath and text field
    import.env$metadataFile$datapath <- ""
    output$metadataText <- renderText({
      paste("No file selected")
    })
    
    # resetting image directory path and text field
    import.env$imageDir <- ""
    output$imageDirtext <- renderText({
      paste("No images directory selected!")
    })
    
    # resetting working directory path and text field
    import.env$Dir <- ""
    output$wdText <- renderText({
      paste("No working directory selected!")
    })
    
    # resetting messages section
    textMsgs(" ")
    
  })
  