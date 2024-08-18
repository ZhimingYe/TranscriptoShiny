library(shiny)
library(shinyjs)
library(bslib)
library(ggplot2)
library(DT)



.PASSWORD<-readRDS("../DEGAnalysis/pwd.rds")
ui <- fluidPage(
  useShinyjs(),
  titlePanel("co-expression Analysis based on NMF"),
  sidebarLayout(
    sidebarPanel(
      div(
        style = "text-align: center;",
        tags$footer(paste0("Version 0.1.2 Author: Zhiming Ye"), style = "font-size: 12px; color: grey;")
      ),
      div(
        style = "text-align: center;",
        tags$footer(paste0("@ YAO LAB"), style = "font-size: 12px; color: grey;")
      ),
      passwordInput("passwd", "Password"),
      hr(),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0("Row: Gene; Column: Sample, Sample name MUST BE unique"),
          style = "font-size: 13px; color: red;"
        )
      ),
      fileInput("countMatrix", "Upload Expression Matrix", accept = c("text/csv","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0(
            "Group Info at least contains TWO columns, the first is sample ID, should be same as expression matrix, the second is grouping information."
          ),
          style = "font-size: 13px; color: black;"
        )
      ),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0(
            "(1)Please don't add any other columns. (2)The first ROW in Excel should start with row names(like SampleID and Group), NOT starts with sample"
          ),
          style = "font-size: 13px; color: red;"
        )
      ),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0("If there is any batch information, please add in the third column"),
          style = "font-size: 13px; color: orange;"
        )
      ),
      fileInput("colData", "Upload Group Information", accept = c("text/csv","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      actionButton("UploadFin", "Access Programme"),
      hr(),
      radioButtons(
        "SEQtypeDfTypeloggedInput",
        "Matrix Type",
        choices = list(
          "TPM (RAW, Suggested)" = "TPMRAW",
          "Counts (RAW, NOT Suggested)" = "COUNTRAW",
          "TPM log-transformed (Suggested)" = "TPMLOG",
          "Not from Sequencing (RAW)" = "LSRAW",
          "Not from Sequencing log-transformed" = "LSLOG"
        ),
        selected = "TPMLOG"
      ),
      sliderInput(
        "NumFilterInput",
        "Genes expressed in N sample were kept",
        min = 1,
        max = 20,
        value = 2
      ),
      radioButtons(
        "SpecInput",
        "Organism",
        choices = list(
          "Human" = "Human",
          "Mouse" = "Mouse",
          "Others" = "Others"
        ),
        selected = "Human"
      ),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0("STRONGly Suggested keep Protein coding genes only"),
          style = "font-size: 13px; color: orange;"
        )
      ),
      conditionalPanel(
        condition = "input.SpecInput != 'Others'",
        radioButtons(
          "FilterPCInput",
          "Filter Gene Type",
          choices = list(
            "Protein coding gene only" = "PC",
            "Not filtering" = "NOT"
          ),
          selected = "PC"
        )),
      hr(),
      checkboxInput(
        "FilterBatchInput",
        "REMOVE batch effect (Might be not correct, only suggest in analysing public data)",
        value = F
      ),
      hr(),
      sliderInput(
        "madselection",
        "Keep N Genes based on MAD",
        min = 2000,
        max = 15000,
        value = 7500
      ),
      sliderInput(
        "testkrange",
        "Test K range (+/-3)",
        min = 4,
        max = 17,
        value = 5
      ),
      actionButton("runDESeq", "Select best cluster nums."),
      hr(),
      sliderInput(
        "TargetKK",
        "Select the elbow K",
        min = 2,
        max = 20,
        value = 4
      ),
      actionButton("runnmf", "Run analysis"),
      hr(),
      downloadButton("downloadData", "Download Gene Contribution"),
      downloadButton("downloadData2", "Download Module-Sample Contribution")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("OutPut", plotOutput("PCAplot",height = "600px")),
        nav_panel("sessionInfo", DTOutput("SesionInfo")),
        full_screen = T,
        wrapper = card_body(height = "1200px")
      )
    )
  )
)
Flag <- F
DEGtable<-data.frame()
source("../DEGAnalysis/TranscriptoShinyLib.R")
server <- function(input, output, session) {
  shinyjs::disable("runDESeq")
  shinyjs::disable("runnmf")
  loadPackages<-function(){
    setProgress(0.35)
    library(dplyr)
    library(sva)
    library(tibble)
    library(readr)
    library(limma)
    library(edgeR)
    library(plyr)
    library(ComplexHeatmap)
    library(readr)
    library(readxl)
    library(RcppML)
    options(RcppML.threads = 1)
  }
  observe({
    if (Flag) {
      shinyjs::enable("downloadData")
      shinyjs::enable("downloadData2")
    } else {
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
    }
  })
  observeEvent(input$UploadFin,{tryCatch({
    if(input$passwd!=.PASSWORD){
      stop("Wrong Password!")
    }
    else{
      showNotification("Welcome!",type = "message")
      loadPackages()
      output$SesionInfo <- renderDT({
        (xfun::session_info() %>% as_tibble())[-c(1,3),]
        
      }, options = list(pageLength = 50, scrollX = T))
    }
    shinyjs::enable("runDESeq")
    shinyjs::disable("runnmf")
    shinyjs::disable("downloadData")
    shinyjs::disable("downloadData2")
    showNotification("Pre-check finish, please continue analysis",type = "message")
  }, error = function(e) {
    showNotification(paste("Error:", e$message), type = "error")
    Flag <- F
  })})
  TargetDF<-NULL
  observeEvent(input$runDESeq, {
    req(input$countMatrix, input$colData)
    
    tryCatch({
      if(input$passwd!=.PASSWORD){
        stop("Wrong Password!")
      }
      shinyjs::disable("runDESeq")
      shinyjs::disable("runnmf")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      countData <- read_file.(input$countMatrix$datapath)
      colData <- read_file.(input$colData$datapath)
      
      if(input$NumFilterInput>ncol(countData)*0.75){
        stop("Not enough column to filter")
      }
      resOrdered<-PrefilterDF(ExprTable = countData,Group = colData,SEQtypeDfTypeloggedInput = input$SEQtypeDfTypeloggedInput,doBatchremove = input$FilterBatchInput,NumFilter = input$NumFilterInput,FilterPC = input$FilterPCInput,Spec = input$SpecInput)
      if(input$SEQtypeDfTypeloggedInput=="COUNTRAW"){
        resOrdered<-edgeR::cpm(resOrdered)
      }
      setProgress(0.5)
      if(input$madselection>nrow(resOrdered)-50){
        stop("Too many gene to keep, please reduce!")
      }
      resOrdered<-CalcMad(resOrdered,input$madselection)
      showNotification("Waiting...This step is very slow...")
      showNotification("It might use 30 min to 2 min")
      showNotification("Please be PATIENT!",type = "message")
      set.seed(2024)
      cv_data <- RcppML::crossValidate(resOrdered, k = (input$testkrange-3):(input$testkrange+3), reps = 7)
      TargetDF<<-resOrdered
      output$PCAplot <- renderPlot({plot(cv_data)+ylab("Mean Squared Error")})
      shinyjs::disable("runDESeq")
      shinyjs::enable("runnmf")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      shinyjs::enable("runDESeq")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      shinyjs::disable("runnmf")
      Flag <<- F
    }
    )%>%withProgress(message = 'Processing...')
  })
  
  observeEvent(input$runnmf, {
    
    tryCatch({
      if(input$passwd!=.PASSWORD){
        stop("Wrong Password!")
      }
      shinyjs::disable("runDESeq")
      shinyjs::disable("runnmf")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      
      MODEL <- RcppML::nmf(TargetDF, k =input$TargetKK, seed=2024)
      output$downloadData <- downloadHandler(
        filename = function() {
          "GeneContribution.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(MODEL@w), file)
        }
      )
      output$downloadData2 <- downloadHandler(
        filename = function() {
          "Contribution2Sample.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(MODEL@h), file)
        }
      )
      shinyjs::enable("runDESeq")
      shinyjs::disable("runnmf")
      shinyjs::enable("downloadData")
      shinyjs::enable("downloadData2")
      # TargetDF<<-resOrdered
      output$PCAplot <- renderPlot({Heatmap((MODEL@h)%>%as.matrix(),name="NMF")})
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      shinyjs::enable("runDESeq")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      shinyjs::disable("runnmf")
      Flag <<- F
    }
    )%>%withProgress(message = 'Processing...')
  })
}
shinyApp(ui = ui, server = server,enableBookmarking="disable")
