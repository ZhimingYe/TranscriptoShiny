library(shiny)
library(shinyjs)
library(bslib)
library(ggplot2)
library(DT)



.PASSWORD<-readRDS("./pwd.rds")
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Differential Expression Analysis"),
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
          "Counts (RAW,Suggested)" = "COUNTRAW",
          "TPM (RAW)" = "TPMRAW",
          "TPM log-transformed" = "TPMLOG",
          "Not from Sequencing (RAW)" = "LSRAW",
          "Not from Sequencing log-transformed" = "LSLOG"
        ),
        selected = "COUNTRAW"
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
            "Protein coding gene only (only supports human and mouse)" = "PC",
            "NCRNA only (only supports human and mouse)" = "NC",
            "Not filtering" = "NOT"
          ),
          selected = "NOT"
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
          "DEA according to Group" = "G",
          "DEA according to Trend" = "T"
        ),
        selected = "G"
      ),
      conditionalPanel(
        condition = "input.ANLtype == 'G'",
        selectInput("selectA", "If according to Group, Select Comparision", choices = NULL)
      ),
      conditionalPanel(
        condition = "input.ANLtype == 'T'",
        selectInput("selectB", "If according to Trend, Select Order", choices = NULL),
        sliderInput(
          "MUZZN",
          "If according to Trend, Select Cluster numbers",
          min = 6,
          max = 15,
          value = 6
        ),
        sliderInput(
          "madselection",
          "Keep N Genes for clustering based on MAD",
          min = 100,
          max = 20000,
          value = 10000
        )
      ),
      hr(),
      conditionalPanel(
        condition = "input.SpecInput != 'Others'",
        checkboxInput("Adgeneann", "Add gene annotation", value = F)),
      actionButton("runDESeq", "Run DEG analysis"),
      downloadButton("downloadData", "Download Differential Genes CSV"),
      conditionalPanel(
        condition = "input.FilterBatchInput | input.SEQtypeDfTypeloggedInput == 'LSRAW'| input.SEQtypeDfTypeloggedInput == 'LSLOG'",
        downloadButton("downloadData2", "Download Batch-Eff Corrected CSV")
      )
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("PCA", plotOutput("PCAplot")),
        nav_panel("Pearson Correlation", plotOutput("PCCCR",height = "600px")),
        nav_panel("DEG analysis", plotOutput("volcanoPlot")),
        nav_panel("DEG table", DTOutput("table")),
        nav_panel("sessionInfo", DTOutput("SesionInfo")),
        full_screen = T,
        wrapper = card_body(height = "1200px")
      )
    )
  )
)
Flag <- F
DEGtable<-data.frame()
source("./TranscriptoShinyLib.R")
server <- function(input, output, session) {
  shinyjs::disable("runDESeq")
  
  loadPackages<-function(){
    # setProgress(0.35)
    library(dplyr)
    library(sva)
    library(tibble)
    library(DESeq2)
    library(readr)
    library(limma)
    library(edgeR)
    library(plyr)
    library(scales)
    library(grid)
    library(FactoMineR)
    library(factoextra)
    library(ComplexHeatmap)
    library(Mfuzz)
    library(RColorBrewer)
    library(readxl)
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
    
    colData <- read_file.(input$colData$datapath)
    GRPINFO<<-GetGroups(colData[,2,drop=T])
    NAME1List<-GRPINFO[["A"]]$GroupOrder
    names(NAME1List)<-GRPINFO[["A"]]$Cpr
    NAME2List<-GRPINFO[["B"]]$GroupOrder
    names(NAME2List)<-GRPINFO[["B"]]$Cpr
    setProgress(0.5)
    updateSelectInput(session, "selectA", choices = NAME1List)
    updateSelectInput(session, "selectB", choices = NAME2List)
    shinyjs::enable("runDESeq")
    shinyjs::disable("downloadData")
    shinyjs::disable("downloadData2")
    setProgress(1)
    showNotification("Pre-check finish, please continue analysis",type = "message")
  }, error = function(e) {
    showNotification(paste("Error:", e$message), type = "error")
    Flag <- F
  })%>%withProgress(message = "Loading...")})
  observeEvent(input$runDESeq, {
    req(input$countMatrix, input$colData)
    
    
    
    tryCatch({
      if(input$passwd!=.PASSWORD){
        stop("Wrong Password!")
      }
      countData <- read_file.(input$countMatrix$datapath)
      colData <- read_file.(input$colData$datapath)
      shinyjs::disable("runDESeq")
      shinyjs::disable("downloadData")
      if(input$NumFilterInput>ncol(countData)*0.75){
        stop("Not enough column to filter")
      }
      resOrdered<-PrefilterDF(ExprTable = countData,Group = colData,SEQtypeDfTypeloggedInput = input$SEQtypeDfTypeloggedInput,doBatchremove = input$FilterBatchInput,NumFilter = input$NumFilterInput,FilterPC = input$FilterPCInput,Spec = input$SpecInput)
      PCAplot<-PlotPCA(Mat = resOrdered,group = colData[,2,drop=T],counts = input$SEQtypeDfTypeloggedInput=="COUNTRAW",loged = T,Ellipse = T,NumOfGenes = 5000)
      CRRplot<-PlotCorr(Mat = resOrdered,group = colData[,2,drop=T],counts = input$SEQtypeDfTypeloggedInput=="COUNTRAW",loged = T,Ellipse = T,NumOfGenes = 5000)
      
      if(input$ANLtype=="T"){
        if(input$SEQtypeDfTypeloggedInput=="COUNTRAW"){
          resOrdered<-edgeR::cpm(resOrdered)
        }
        TgtFactor<-factor(colData[,2,drop=T],levels = strsplit(GRPINFO[["B"]]$Cpr[GRPINFO[["B"]]$GroupOrder==input$selectB], "-")[[1]])
        
        MfuzzMat<-apply(resOrdered,1,function(x)tapply(x,TgtFactor,mean))
        MfuzzMat<-MfuzzMat%>%t()%>%as.data.frame()
        setProgress(0.5)
        showNotification("Preparing for Mfuzz Clustering...")
        MfuzzMat<-CalcMad(MfuzzMat,input$madselection)
        library(SummarizedExperiment)
        MfuzzDs1<-new('ExpressionSet',exprs=MfuzzMat%>%as.matrix())
        MfuzzDs1 <- standardise(MfuzzDs1)
        set.seed(2024)
        showNotification("Mfuzz Clustering...")
        cl <<- mfuzz(MfuzzDs1,c=input$MUZZN,m= mestimate(MfuzzDs1))
        DEGtable<-cl[["cluster"]]%>%as.data.frame()
        DEGtable<-DEGtable%>%rownames_to_column(var="GID")
        colnames(DEGtable)[2]<-"Belong to Cluster"
        if(input$Adgeneann){
          DEGtable<-Get0ID(DEGtable,SpecInput = input$SpecInput)
        }
        DEGtable<<-DEGtable
        MfuzzDs1a<<-MfuzzDs1
        output$volcanoPlot <- renderPlot({mfuzz.plot(MfuzzDs1a, cl = cl,time.labels =strsplit(GRPINFO[["B"]]$Cpr[GRPINFO[["B"]]$GroupOrder==input$selectB], "-")[[1]],  mfrow = c(3, 5), new.window = FALSE, colo = rev(brewer.pal(11, "RdBu")))})
        Flag <<- T
      }
      if(input$ANLtype=="G"){
        if(input$SEQtypeDfTypeloggedInput%in%c("COUNTRAW")){
          setProgress(0.5)
          DEGtable<-GenDEGtable(Expr = resOrdered,GroupInput = colData[,2,drop=T],CprString = GRPINFO[["A"]]$Cpr[GRPINFO[["A"]]$GroupOrder==input$selectA])
          if(input$Adgeneann){
            DEGtable<-Get0ID(DEGtable,SpecInput = input$SpecInput)
          }
          DEGtable<<-DEGtable
        }
        else{
          setProgress(0.5)
          DEGtable<-GenDEGtable.Norm(Expr = resOrdered,GroupInput = colData[,2,drop=T],CprString = GRPINFO[["A"]]$Cpr[GRPINFO[["A"]]$GroupOrder==input$selectA])
          DEGtable<<-DEGtable
        }
        dfa<-DEGtable%>%dplyr::filter(Package=="limma")
        dfa$significant <- with(dfa, ifelse(Padj < 0.05 & abs(logFC) > 1, "Significant", ifelse(PValue < 0.05 & abs(logFC) > 0.25,"Might be significant","Not Significant")))
        px <- ggplot(dfa, aes(x=logFC, y=-log10(Padj))) +
          geom_point(aes(color=significant), alpha=0.8, size=1.75) +
          scale_color_manual(values=c("#5ecd72", "#88a8d3","#665571")) +
          theme_minimal() +
          labs(title="Volcano Plot",
               x="Log2 Fold Change",
               y="-Log10 adj. p-value") +
          theme(legend.position = "top")+theme_bw()
        showNotification("Analysis Finished!",type = "message")
        output$volcanoPlot <- renderPlot({px})
        Flag <<- T
      }
      output$PCAplot <- renderPlot({
        if (Flag) {
          PCAplot
        }
        
      })
      output$PCCCR <- renderPlot({
        if (Flag) {
          CRRplot
        }})
      output$table <- renderDT({
        if (Flag) {
          datatable(as.data.frame(DEGtable), options =  list(pageLength = 50, scrollX = T))
        }
      })
      output$downloadData <- downloadHandler(
        filename = function() {
          "differential_genes.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(DEGtable), file)
        }
      )
      output$downloadData2 <- downloadHandler(
        filename = function() {
          "ExprMatCorr.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(resOrdered), file)
        }
      )
      shinyjs::enable("downloadData")
      shinyjs::enable("runDESeq")
      if(input$FilterBatchInput|input$SEQtypeDfTypeloggedInput%in%c("LSRAW","LSLOG")){
        shinyjs::enable("downloadData2")
      }
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      shinyjs::enable("runDESeq")
      Flag <<- F
    }
    )%>%withProgress(message = 'Processing...')
    
  })
  
  
}


shinyApp(ui = ui, server = server,enableBookmarking="disable")
