library(shiny)
library(shinyjs)
library(bslib)
library(ggplot2)
library(DT)




.PASSWORD<-readRDS("../DEGAnalysis/pwd.rds")
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Revealing biological functions by Sample Scoring"),
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
      fileInput("GMTFile", "Provide Self-built gene set (GMT files, Optional)", accept = c(".gmt")),
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
      conditionalPanel(
        condition = "input.SpecInput != 'Others'",
        radioButtons(
          "FilterPCInput",
          "Filter Gene Type",
          choices = list(
            "Protein coding gene only" = "PC",
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
          "PROGENy model pathway activity" = "P_DecoupleR",
          "CollecTRI TF Infering" = "C_DecoupleR",
          "xCell Immune Cell" = "Imm_xCell",
          "tumor HALLMARK"="Hallmark_ssGSEA",
          "KEGG_LEGACY ssGSEA"="KeggOld_ssGSEA",
          "KEGG_MEDICUS ssGSEA"="KeggNew_ssGSEA",
          "BioCarta ssGSEA"="BioCarta_ssGSEA",
          "WikiPathways ssGSEA"="WikiPathways_ssGSEA",
          "PID ssGSEA"="PID_ssGSEA",
          "oncogenic sig ssGSEA"="OncoSig_ssGSEA",
          "self-built gene set ssGSEA"="SB_ssGSEA",
          "self-built gene set ORA"="SB_ORA"
        ),
        selected = "P_DecoupleR"
      ),
      actionButton("runDESeq", "Run analysis"),
      downloadButton("downloadData", "Download Result")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("Result Preview", plotOutput("PCAplot",height = "1150px")),
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
EnsemblConvertor<-function(x,org){
  tpmset<-x
  if(org=="Human"){
    TPMgeneinfo<-bitr(rownames(tpmset),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
  }
  else{
    TPMgeneinfo<-bitr(rownames(tpmset),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
  }
  tpmsetNAME<-tpmset%>%as.data.frame()%>%rownames_to_column(var="ENSEMBL")
  tpmsetNAME<-tpmsetNAME%>%left_join(TPMgeneinfo)
  tpmsetNAME<-tpmsetNAME[!is.na(tpmsetNAME$SYMBOL),]
  tpmsetNAMESUM<-aggregate(tpmsetNAME[,-which(colnames(tpmsetNAME)%in%c("ENSEMBL","SYMBOL"))],by=list(SYMBOL=tpmsetNAME$SYMBOL),FUN="sum")
  tpmsetNAMESUM<-tpmsetNAMESUM%>%as.data.frame()%>%column_to_rownames(var="SYMBOL")
  return(tpmsetNAMESUM)
}

server <- function(input, output, session) {
  shinyjs::disable("runDESeq")
  
  
  loadPackages<-function(){
    setProgress(0.35)
    library(dplyr)
    library(sva)
    library(tibble)
    library(readr)
    library(edgeR)
    library(plyr)
    library(scales)
    library(grid)
    library(ComplexHeatmap)
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
      if(input$IDType=="ENS"){
        resOrdered<-EnsemblConvertor(resOrdered,org = input$SpecInput)
      }
      needSTD<-F
      load("../DEGAnalysis/AnnoPwDBsV2.RData")
      if(input$ANLtype=="P_DecoupleR"){
        if(input$SpecInput=="Human"){
          net<-progenyNet.human
          rownames(resOrdered)<-stringr::str_to_upper(rownames(resOrdered))
        }
        else{
          net<-progenyNet.mouse
          rownames(resOrdered)<-stringr::str_to_title(rownames(resOrdered))
        }
        setProgress(0.5)
        sample_acts <- run_mlm(mat=resOrdered, net=net, .source='source', .target='target',.mor='weight', minsize = 5)
        setProgress(0.75)
        library(reshape2)
        sample_acts02<-dcast(source~condition,value.var = "score",data = sample_acts)
        gsva_es<-sample_acts02%>%column_to_rownames(var="source")
      }
      if(input$ANLtype=="C_DecoupleR"){
        if(input$SpecInput=="Human"){
          net<-collectriNet.human
          rownames(resOrdered)<-stringr::str_to_upper(rownames(resOrdered))
        }
        else{
          net<-collectriNet.mouse
          rownames(resOrdered)<-stringr::str_to_title(rownames(resOrdered))
        }
        setProgress(0.5)
        sample_acts <- run_ulm(mat=resOrdered, net=net, .source='source', .target='target',.mor='weight', minsize = 5)
        setProgress(0.75)
        library(reshape2)
        sample_acts02<-dcast(source~condition,value.var = "score",data = sample_acts)
        gsva_es<-sample_acts02%>%column_to_rownames(var="source")
      }
      
      if(grepl("ssGSEA",input$ANLtype)|grepl("xCell",input$ANLtype)){
        setProgress(0.5)
        rownames(resOrdered)<- stringr::str_to_upper(rownames(resOrdered))
        needSTD<-T
      }
      if(grepl("ssGSEA",input$ANLtype)){
        setProgress(0.5)
        resOrdered<-as.matrix(resOrdered)
      }
      if(input$ANLtype=="Imm_xCell"){
        setProgress(0.5)
        gsva_es<-xCellAnalysis(resOrdered)
      }
      if(input$ANLtype=="Hallmark_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, h.all)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="KeggOld_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c2.cp.kegg_legacy)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="KeggNew_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c2.cp.kegg_medicus)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="BioCarta_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c2.cp.biocarta)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="WikiPathways_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c2.cp.wikipathways)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="PID_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c2.cp.pid)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="OncoSig_ssGSEA"){
        setProgress(0.5)
        sp1 <- ssgseaParam(resOrdered, c6.all)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="SB_ssGSEA"){
        setProgress(0.5)
        GMTobj <- getGmt(input$GMTFile$datapath)
        if(length(GMTobj)>850){
          stop("Too large GMT file, please use your self computer!")
        }
        sp1 <- ssgseaParam(resOrdered, GMTobj)
        gsva_es <- gsva(sp1)
      }
      if(input$ANLtype=="SB_ORA"){
        setProgress(0.5)
        GMTobj <- read.gmt(input$GMTFile$datapath)
        
        if(length(names(table(GMTobj$term)))>=100){
          stop("Too large GMT file, please use your self computer!")
        }
        GMTobj<-as.data.frame(GMTobj)
        # [1] "term" "gene"
        GMTobj$GT<-paste0(GMTobj$term,GMTobj$gene)
        GMTobj<-GMTobj[!duplicated(GMTobj$GT),]
        colnames(GMTobj)[c(1,2)]<-c("source","target")
        sample_acts <- run_ora(mat=resOrdered, net=GMTobj, .source='source', .target='target', minsize = 5,n_up=round(0.15*nrow(resOrdered)))
        showNotification("If result is odd, please use filtering protein-coding gene function or use ssGSEA instead. ")
        setProgress(0.75)
        library(reshape2)
        sample_acts02<-dcast(source~condition,value.var = "score",data = sample_acts)
        gsva_es<-sample_acts02%>%column_to_rownames(var="source")
        # sp1 <- ssgseaParam(resOrdered, GMTFile)
        # gsva_es <- gsva(sp1)
      }
      output$PCAplot <- renderPlot({
        if(needSTD){
          tryCatch({Heatmap(std_(gsva_es),name="ES")}, error = function(e) {Heatmap((gsva_es),name="ES")})
        }
        else{
          Heatmap((gsva_es),name="ES")
        }
          
        
      })
      output$downloadData <- downloadHandler(
        filename = function() {
          "ES.csv"
        },
        content = function(file) {
          write.csv(as.data.frame(gsva_es), file)
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
