library(shiny)
library(shinyjs)
library(dplyr)
library(sva)
library(tibble)
library(bslib)
# library(DESeq2)
library(ggplot2)
library(DT)
library(readr)
# library(limma)
# library(DESeq2)
library(ggplot2)
library(edgeR)
# library(limma)
library(plyr)
library(scales)
library(grid)
# library(FactoMineR)
# library(factoextra)
#gtools
library(ComplexHeatmap)
# library(Mfuzz)
library(RColorBrewer)
library(readr)
library(readxl)
library(xCell)
library(decoupleR)
library(reshape2)

library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(gson)
library(ReactomePA)
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Enrichment analysis"),
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
        tags$footer(paste0("For ORA: The first column name is 'Gene'. If there is group information provide them as the second column. Name as 'Group'. "), style = "font-size: 13px; color: blue;")
      ),
      div(
        style = "text-align: left;",
        tags$footer(paste0("For GSEA: The first column name is 'Gene'. The second is 'Weight'"), style = "font-size: 13px; color: red;")
      ),
      fileInput("countMatrix", "Upload Gene List", accept = c("text/csv","application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")),
      fileInput("GMTFile", "Provide Self-built gene set (GMT files, Optional)", accept = c(".gmt")),
      actionButton("UploadFin", "Access Programme"),
      hr(),
      div(
        style = "text-align: left;",
        tags$footer(paste0("WARNING: In enrichment analysis, all analysis is SLOW, please be patient."), style = "font-size: 15px; color: red;")
      ),
      radioButtons(
        "SpecInput",
        "Organism",
        choices = list(
          "Human" = "Human",
          "Mouse" = "Mouse"
        ),
        selected = "Human"
      ),
      radioButtons(
        "IDType",
        "Gene ID Type",
        choices = list(
          "Symbol" = "SYM",
          "Ensembl ID" = "ENS"
        ),
        selected = "SYM"
      ),
      hr(),
      radioButtons(
        "MultiGroup",
        "Analysis Type",
        choices = list(
          "ORA" = "SG",
          "Multi-group ORA" = "MG",
          "GSEA" = "GSEA"
        ),
        selected = "SG"
      ),
      radioButtons(
        "ANLtype",
        "Use Database",
        choices = list(
          "GO Biological Process" = "gobp",
          "GO Cell Compounent" = "gocc",
          "GO Molecular Function" = "gomf",
          "KEGG (typical)" = "kegg",
          "MKEGG" = "mkegg", 
          "Reactome" = "ra",
          "Disease Ontology"="do",
          "tumor HALLMARK"="hal",
          "KEGG_LEGACY via MSIGDB"="kl",
          "KEGG_MEDICUS via MSIGDB"="km",
          "BioCarta"="bc",
          "WikiPathways"="wp",
          "PID"="pid",
          "oncogenic sig (MSIGDB C6)"="os",
          "Cell Type sig (MSIGDB C8)"="cs",
          "self-built gene set"="SB"
        ),
        selected = "ra"
      ),
      actionButton("runDESeq", "Run analysis"),
      downloadButton("downloadData", "Download Result"),
      downloadButton("downloadData2", "Download RDS")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("Result Preview", plotOutput("PCAplot",height = "800px",width = "600px")),
        nav_panel("enrichment DF", DTOutput("CorrOut")),
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
    if(input$passwd!=""){
      stop("Wrong Password!")
    }
    else{
      showNotification("Welcome!",type = "message")
      output$SesionInfo <- renderDT({
        (xfun::session_info() %>% as_tibble())[-c(1,3),]
        
      }, options = list(pageLength = 50, scrollX = T))
    }
    shinyjs::enable("runDESeq")
    shinyjs::disable("downloadData")
    shinyjs::disable("downloadData2")
    showNotification("Pre-check finish, please continue analysis",type = "message")
  }, error = function(e) {
    showNotification(paste("Error:", e$message), type = "error")
    Flag <- F
  })})
  observeEvent(input$runDESeq, {
    req(input$countMatrix)
    
    tryCatch({
      if(input$passwd!=""){
        stop("Wrong Password!")
      }
      shinyjs::disable("runDESeq")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      countData <- read_file.(input$countMatrix$datapath)
      colnames(countData)<-stringr::str_to_title(colnames(countData))
      
      if(input$IDType=="ENS"){
        colnames(countData)[which(colnames(countData)=="Gene")]<-"ENSEMBL"
        countData<-ConvertEnsembl2Symbol(countData,input$SpecInput)
        
      }
      else{
        colnames(countData)[which(colnames(countData)=="Gene")]<-"SYMBOL"
        if(input$SpecInput=="Mouse"){
          countData<-MouseSymbol2Human(countData)
        }
      }
      # countData2<<-countData
      GSt<-NULL
      library(org.Hs.eg.db)
      notmsigdblist<-c("gobp","gocc","gomf","kegg","mkegg","ra","do")
      if(input$MultiGroup=="MG"){
        if(input$ANLtype%in%notmsigdblist){
          GSt<-BuildMultigroupDEGlist(countData$SYMBOL,countData$Group,UniversalDB = F)
        }
        else{
          GSt<-BuildMultigroupDEGlist(countData$SYMBOL,countData$Group,UniversalDB = T)
        }
      }
      if(input$MultiGroup=="GSEA"){
        if(input$ANLtype%in%notmsigdblist){
          GSt<-Ranked.GS(countData$SYMBOL,countData$Weight,UniversalDB = F)
        }
        else{
          GSt<-Ranked.GS(countData$SYMBOL,countData$Weight,UniversalDB = T)
        }
      }
      if(input$MultiGroup=="SG"){
        GSt<-countData$SYMBOL
      }
      load("../DEGAnalysis/AnnoPwDBsV2_4CP.RData")
      if(input$ANLtype=="gobp"){
        setProgress(0.5)
        gsva_es <- doGO(MultiGroup = input$MultiGroup,GS=GSt,ont = "BP",PVal = 0.1,QVal = 0.5)
      }
      if(input$ANLtype=="gocc"){
        setProgress(0.5)
        gsva_es <- doGO(MultiGroup = input$MultiGroup,GS=GSt,ont = "CC",PVal = 0.1,QVal = 0.5)
      }
      if(input$ANLtype=="gomf"){
        setProgress(0.5)
        gsva_es <- doGO(MultiGroup = input$MultiGroup,GS=GSt,ont = "MF",PVal = 0.1,QVal = 0.5)
      }
      if(input$ANLtype=="kegg"){
        setProgress(0.5)
        gsva_es <- doKEGG(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,GSONFILE = KEGGgson)
      }
      if(input$ANLtype=="mkegg"){
        setProgress(0.5)
        gsva_es <- doMKEGG(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,GSONFILE = MKEGGgson)
      }
      if(input$ANLtype=="ra"){
        setProgress(0.5)
        gsva_es <- doRA(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5)
      }
      if(input$ANLtype=="do"){
        setProgress(0.5)
        gsva_es <- doDO(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5)
      }
      if(input$ANLtype=="hal"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = h.all)
      }
      if(input$ANLtype=="kl"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c2.cp.kegg_legacy)
      }
      if(input$ANLtype=="km"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c2.cp.kegg_medicus)
      }
      if(input$ANLtype=="bc"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c2.cp.biocarta)
      }
      if(input$ANLtype=="wp"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c2.cp.wikipathways)
      }
      if(input$ANLtype=="pid"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c2.cp.pid)
      }
      if(input$ANLtype=="os"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c6.all)
      }
      if(input$ANLtype=="cs"){
        setProgress(0.5)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = c8.all)
      }
      if(input$ANLtype=="SB"){
        setProgress(0.5)
        GMTobj <- read.gmt(input$GMTFile$datapath)
        gsva_es <- doUniversal(MultiGroup = input$MultiGroup,GS=GSt,PVal = 0.1,QVal = 0.5,DB = GMTobj)
      }
      
      output$PCAplot <- renderPlot({
          tryCatch({
            if(input$MultiGroup!="GSEA"){
              enrichplot::dotplot(gsva_es)
            }
            else{
              enrichplot::ridgeplot(gsva_es)
            }
          }, error = function(e) {showNotification("ERROR",type = "error")})})
      output$CorrOut <- renderDT({
        gsva_es%>%as.data.frame()
        
      }, options = list(pageLength = 50, scrollX = T))
      output$downloadData <- downloadHandler(
        filename = function() {
          "Enrichment.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(gsva_es), file)
        }
      )
      output$downloadData2 <- downloadHandler(
        filename = function() {
          "Enrichment.rds"
        },
        content = function(file) {
          saveRDS(gsva_es, file)
        }
      )
      shinyjs::enable("runDESeq")
      shinyjs::enable("downloadData")
      shinyjs::enable("downloadData2")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      shinyjs::enable("runDESeq")
      shinyjs::disable("downloadData")
      shinyjs::disable("downloadData2")
      Flag <<- F
    }
    )%>%withProgress(message = 'Processing...')
  })
}
shinyApp(ui = ui, server = server,enableBookmarking="disable")
