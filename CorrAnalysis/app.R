library(shiny)
library(shinyjs)
library(ggplot2)
library(DT)

library(bslib)


.PASSWORD<-readRDS("../DEGAnalysis/pwd.rds")
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Batch Correlation Analysis of Numeric Matrix"),
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
      radioButtons(
        "ANLtype",
        "Analysis Type",
        choices = list(
          "Pearson Correlation" = "PCC",
          "Spearman Correlation" = "SPC"
        ),
        selected = "PCC"
      ),
      textInput("Target","Target Gene"),
      actionButton("runDESeq", "Run analysis"),
      downloadButton("downloadData", "Download Result")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("Result Preview", DTOutput("CorrOut")),
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
  loadPackages<-function(){
    # setProgress(0.35)
    library(tibble)
    library(readr)
    library(limma)
    library(dplyr)
    library(sva)
    library(edgeR)
    library(plyr)
    library(readr)
    library(readxl)
  }
  observe({
    if (Flag) {
      shinyjs::enable("downloadData")
    } else {
      shinyjs::disable("downloadData")
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
    shinyjs::disable("downloadData")
    showNotification("Pre-check finish, please continue analysis",type = "message")
  }, error = function(e) {
    showNotification(paste("Error:", e$message), type = "error")
    Flag <- F
  })})
  observeEvent(input$runDESeq, {
    req(input$countMatrix, input$colData)
    
    tryCatch({
      if(input$passwd!=.PASSWORD){
        stop("Wrong Password!")
      }
      shinyjs::disable("runDESeq")
      shinyjs::disable("downloadData")
      countData <- read_file.(input$countMatrix$datapath)
      colData <- read_file.(input$colData$datapath)
      
      if(input$NumFilterInput>ncol(countData)*0.75){
        stop("Not enough column to filter")
      }
      resOrdered<-PrefilterDF(ExprTable = countData,Group = colData,SEQtypeDfTypeloggedInput = input$SEQtypeDfTypeloggedInput,doBatchremove = input$FilterBatchInput,NumFilter = input$NumFilterInput,FilterPC = input$FilterPCInput,Spec = input$SpecInput)
      if(input$SEQtypeDfTypeloggedInput=="COUNTRAW"){
        resOrdered<-edgeR::cpm(resOrdered)
      }
      
      if(!input$Target%in%rownames(resOrdered)){
        stop("Target not in Expression matrix")
      }
      if(input$ANLtype=="PCC"){
        CorrOutDF<-CorrEstimate0(resOrdered,TargetSet = input$Target,UseMethod = "pearson")
      }
      else{
        CorrOutDF<-CorrEstimate0(resOrdered,TargetSet = input$Target,UseMethod = "spearman")
      }
      output$CorrOut <- renderDT({
        CorrOutDF
        
      }, options = list(pageLength = 50, scrollX = T))
      
      output$downloadData <- downloadHandler(
        filename = function() {
          "Correlation.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(CorrOutDF), file)
        }
      )
      shinyjs::enable("runDESeq")
      shinyjs::enable("downloadData")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      shinyjs::enable("runDESeq")
      Flag <<- F
    }
    )%>%withProgress(message = 'Processing...')
  })
}
shinyApp(ui = ui, server = server,enableBookmarking="disable")
