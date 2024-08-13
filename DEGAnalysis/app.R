library(shiny)
library(shinyjs)
library(bslib)
library(DESeq2)
library(ggplot2)
library(DT)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Differential Expression Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("countMatrix", "Upload Expression Matrix (CSV)", accept = "text/csv"),
      fileInput("colData", "Upload Group Information (CSV)", accept = "text/csv"),
      selectInput("selectA", "Select Comparasion", choices = NULL),
      actionButton("runDESeq", "Run DEG analysis"),
      downloadButton("downloadData", "Download Differential Genes CSV")
    ),
    mainPanel(
      navset_card_underline(
        nav_panel("PCA", plotOutput("PCAplot")),
        nav_panel("Pearson Correlation", plotOutput("PCCCR")),
        nav_panel("DEG analysis",plotOutput("volcanoPlot")),
        nav_panel("DEG table",DTOutput("table"))
      )
    )
  )
)
Flag<-F

server <- function(input, output, session) {
  observe({
    if (Flag) {
      shinyjs::enable("downloadData")
    } else {
      shinyjs::disable("downloadData")
    }
  })
  colData<-data.frame(Group=NULL)
  observe({
    if (!is.null(colData)) {
      updateSelectInput(session, "selectA", choices = unique(colData$Group))
    }
  })
  observeEvent(input$runDESeq, {
    req(input$countMatrix, input$colData)
    tryCatch({
      countData <- read.csv(input$countMatrix$datapath, row.names = 1)
      colData <- read.csv(input$colData$datapath, row.names = 1)
      
      if(!"condition" %in% colnames(colData)) {
        stop("The group information file must contain a 'condition' column.")
      }
      
      dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
      dds <- DESeq(dds)
      res <- results(dds)
      resOrdered <- res[order(res$padj), ]
      
      Flag<-T
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      Flag<-F
    })
    
  })
  
  
  output$volcanoPlot <- renderPlot({
    if(Flag){
      ggplot(resOrdered, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(alpha = 0.4) +
        xlab("Log2 Fold Change") +
        ylab("-Log10 Adjusted P-value") +
        theme_minimal()
    }
    
  })
  
  output$table <- renderDT({
    if(Flag){
      datatable(as.data.frame(resOrdered), options = list(pageLength = 10))
    }
    
  })
  
  output$downloadData <- downloadHandler(
      filename = function() { "differential_genes.csv" },
      content = function(file) {
        write.csv(as.data.frame(resOrdered), file)
      }
  )
  
  
  
}

shinyApp(ui = ui, server = server)