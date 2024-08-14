library(shiny)
library(shinyjs)
library(tibble)
library(bslib)
library(DESeq2)
library(ggplot2)
library(DT)

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
      radioButtons(
        "radio",
        "Matrix Type",
        choices = list(
          "Counts" = 1,
          "TPM (RAW)" = 2,
          "TPM log-transformed" = 3,
          "Not from Sequencing (RAW)" = 4,
          "Not from Sequencing log-transformed" = 5
        ),
        selected = 1
      ),
      sliderInput(
        "sliderRows",
        "Genes expressed in N sample were keet",
        min = 1,
        max = 1,
        value = 1
      ),
      radioButtons(
        "radio3",
        "Filter Gene Type",
        choices = list(
          "Protein coding gene only (only supports human and mouse)" = "G",
          "NCRNA only (only supports human and mouse)" = "T",
          "Not filtering" = "N"
        ),
        selected = "N"
      ),
      hr(),
      checkboxInput(
        "checkbox",
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
      checkboxInput("checkbox", "Add gene annotation", value = TRUE),
      div(
        style = "text-align: center;",
        tags$text(paste0("Annotation only supports human and mouse."), style = "font-size: 12px; color: black;")
      ),
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



PrefilterTable<-function(df1,df2,mattype,ngenes,filtgene,batrm){
  qzhFPKM[,1,drop=T]%>%typeof()
  sum(duplicated(qzhFPKM[,1,drop=T]))
}



server <- function(input, output, session) {
  observe({
    if (Flag) {
      shinyjs::enable("downloadData")
    } else {
      shinyjs::disable("downloadData")
    }
  })
  colData <- data.frame(Group = NULL)
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
      
      if (!"condition" %in% colnames(colData)) {
        stop("The group information file must contain a 'condition' column.")
      }
      
      dds <- DESeqDataSetFromMatrix(
        countData = countData,
        colData = colData,
        design = ~ condition
      )
      dds <- DESeq(dds)
      res <- results(dds)
      resOrdered <- res[order(res$padj), ]
      
      Flag <- T
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      Flag <- F
    })
    
  })
  
  
  output$volcanoPlot <- renderPlot({
    if (Flag) {
      ggplot(resOrdered, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(alpha = 0.4) +
        xlab("Log2 Fold Change") +
        ylab("-Log10 Adjusted P-value") +
        theme_minimal()
    }
    
  })
  
  output$table <- renderDT({
    if (Flag) {
      datatable(as.data.frame(resOrdered), options = list(pageLength = 10))
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

shinyApp(ui = ui, server = server)