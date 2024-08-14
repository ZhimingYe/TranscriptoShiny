library(shiny)
library(shinyjs)
library(dplyr)
library(sva)
library(tibble)
library(bslib)
library(DESeq2)
library(ggplot2)
library(DT)
library(readr)
library(limma)
library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(FactoMineR)
library(factoextra)
library(ggbiplot)
library(edgeR)
library(ComplexHeatmap)
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Differential Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      div(
        style = "text-align: center;",
        tags$footer(paste0("Version 1.0 Author: Zhiming Ye"), style = "font-size: 12px; color: grey;")
      ),
      div(
        style = "text-align: center;",
        tags$footer(paste0("@ YAO LAB"), style = "font-size: 12px; color: grey;")
      ),
      passwordInput("passwd", "Password"),
      actionButton("loginBtn", "Login"),
      hr(),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0("Please save as CSV using Excel (File->save as comma seperated)"),
          style = "font-size: 13px; color: black;"
        )
      ),
      div(
        style = "text-align: left;",
        tags$footer(paste0("WARNING: DONOT use UTF8 csv!"), style = "font-size: 13px; color: red;")
      ),
      div(
        style = "text-align: left;",
        tags$footer(
          paste0("Row: Gene; Column: Sample, Sample name MUST BE unique"),
          style = "font-size: 13px; color: red;"
        )
      ),
      fileInput("countMatrix", "Upload Expression Matrix (CSV)", accept = "text/csv"),
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
      fileInput("colData", "Upload Group Information (CSV)", accept = "text/csv"),
      actionButton("UploadFin", "Uploaded Confirm"),
      radioButtons(
        "SEQtypeDfTypeloggedInput",
        "Matrix Type",
        choices = list(
          "Counts" = "COUNTRAW",
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
        "radio2",
        "Analysis Type",
        choices = list(
          "DEA according to Group" = "G",
          "DEA according to Trend" = "T"
        ),
        selected = "G"
      ),
      conditionalPanel(
        condition = "input.radio2 == 'G'",
        selectInput("selectA", "If according to Group, Select Comparision", choices = NULL)
      ),
      conditionalPanel(
        condition = "input.radio2 == 'T'",
        selectInput("selectB", "If according to Trend, Select Order", choices = NULL),
        sliderInput(
          "sliderRows",
          "If according to Trend, Select Cluster numbers",
          min = 1,
          max = 10,
          value = 1
        ),
        checkboxInput(
          "checkbox",
          "If according to Trend, Perform prefiltering based on DEA?",
          value = TRUE
        ),
        checkboxInput("checkbox", "If according to Trend, adjust P?", value = TRUE)
      ),
      hr(),
      conditionalPanel(
        condition = "input.SpecInput != 'Others'",
        checkboxInput("checkbox", "Add gene annotation", value = F)),
      actionButton("runDESeq", "Run DEG analysis"),
      downloadButton("downloadData", "Download Differential Genes CSV")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("PCA", plotOutput("PCAplot")),
        nav_panel("Pearson Correlation", plotOutput("PCCCR")),
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
source("~/Documents/4Fun/TranscriptoShiny/TranscriptoShinyLib.R")


# PrefilterTable<-function(df1,df2,mattype,ngenes,filtgene,batrm){
#   qzhFPKM[,1,drop=T]%>%typeof()
#   sum(duplicated(qzhFPKM[,1,drop=T]))
# }
MaxCol<-function(ExprTable){
  if(ncol(ExprTable)>=4){
    return(round(ncol(ExprTable)/1.25))
  }
  else{
    stop("Not enough Samples!")
  }
}


server <- function(input, output, session) {
  shinyjs::disable("runDESeq")
  observe({
    if (Flag) {
      shinyjs::enable("downloadData")
    } else {
      shinyjs::disable("downloadData")
    }
  })
  observeEvent(input$UploadFin,{tryCatch({
    colData <- read_csv(input$colData$datapath)
    GRPINFO<-GetGroups(colData[,2,drop=T])
    NAME1List<-GRPINFO[["A"]]$GroupOrder
    names(NAME1List)<-GRPINFO[["A"]]$Cpr
    NAME2List<-GRPINFO[["B"]]$GroupOrder
    names(NAME2List)<-GRPINFO[["B"]]$Cpr
    updateSelectInput(session, "selectA", choices = NAME1List)
    updateSelectInput(session, "selectB", choices = NAME2List)
    shinyjs::enable("runDESeq")
  }, error = function(e) {
    showNotification(paste("Error:", e$message), type = "error")
    Flag <- F
  })})
  observeEvent(input$runDESeq, {
    req(input$countMatrix, input$colData)
    
    
    
    tryCatch({
      countData <- read_csv(input$countMatrix$datapath)
      colData <- read_csv(input$colData$datapath)
      if(input$NumFilterInput>ncol(countData)*0.75){
        stop("Not enough column to filter")
      }
      resOrdered<-PrefilterDF(ExprTable = countData,Group = colData,SEQtypeDfTypeloggedInput = input$SEQtypeDfTypeloggedInput,doBatchremove = input$FilterBatchInput,NumFilter = input$NumFilterInput,FilterPC = input$FilterPCInput,Spec = input$SpecInput)
      PCAplot<-PlotPCA(Mat = resOrdered,group = colData[,2,drop=T],counts = input$SEQtypeDfTypeloggedInput=="COUNTRAW",loged = T,Ellipse = T,NumOfGenes = 5000)
      CRRplot<-PlotCorr(Mat = resOrdered,group = colData[,2,drop=T],counts = input$SEQtypeDfTypeloggedInput=="COUNTRAW",loged = T,Ellipse = T,NumOfGenes = 5000)
      
      Flag <- T
      shinyjs::enable("downloadData")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      Flag <- F
    })
    output$PCAplot <- renderPlot({
      if (Flag) {
        PCAplot
      }
      
    })
    output$PCCCR <- renderPlot({
      if (Flag) {
        CRRplot
      }
      
    })
    output$table <- renderDT({
      if (Flag) {
        datatable(as.data.frame(resOrdered), options =  list(pageLength = 50, scrollX = T))
      }
      
    })
  })
  
  
  output$volcanoPlot <- renderPlot({
    if (Flag) {
      # ggplot(resOrdered, aes(x = log2FoldChange, y = -log10(padj))) +
      #   geom_point(alpha = 0.4) +
      #   xlab("Log2 Fold Change") +
      #   ylab("-Log10 Adjusted P-value") +
      #   theme_minimal()
    }
    
  })
  
  
  output$SesionInfo <- renderDT({
    (xfun::session_info() %>% as_tibble())[-c(1:3), ]
  }, options = list(pageLength = 50, scrollX = T))
  output$downloadData <- downloadHandler(
    filename = function() {
      "differential_genes.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(resOrdered), file)
    }
  )
}

# df<-read_csv ("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/工作簿2.csv")[,2,drop=T]
# 
# PCAplot<-PlotPCA(Mat = resOrdered,group = df,counts = T,loged = T,Ellipse = T,NumOfGenes = 5000)


# 
# sum(((read_csv ("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/工作簿2.csv")[,2,drop=T])%>%table()%>%as.numeric())>=3)>=2
# 
# DFA<-read_csv("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/GSE223056_geo_counts_副本.csv")
# duplicated(DFA$gene)%>%sum()
# 
# read_csv ("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/工作簿2.csv")[, 1, drop = T]
# resOrdered<-PrefilterDF(ExprTable = read_csv("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/GSE223056_geo_counts_副本.csv"),Group =read_csv ("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/工作簿2.csv"),SEQtypeDfTypeloggedInput = "TPMRAW",doBatchremove = T,NumFilter = 2,FilterPC = "PC",Spec ="Mouse")

# 启动一个包含分组信息的vector

shinyApp(ui = ui, server = server)
# 
# 
# 
# library(gtools)
# 
# # 定义分组信息
# group <- c("A", "B")
# 
# # 生成所有可能的排列
# all_permutations <- permutations(n = length(group), r = length(group), v = group)
# 
# all_permutations <- apply(all_permutations, 1, paste, collapse = "")
# 
# dd<-list(`A`=data.frame(group=c("A", "B")))
# dd$A


colData <- (read_csv ("~/Documents/A_onGoing/MemGradFillin2Days/YMJ_LK99/QZH_Glia_CELLs_inLUNG/工作簿2.csv"))
GRPINFO<-GetGroups(colData[,2,drop=T])

un