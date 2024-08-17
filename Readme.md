# Transcriptomic analysis on Shiny

# Purpose

I come from a wet lab where most members are not familiar with using R programming for sequencing data analysis. Therefore, this [Shiny](https://shiny.posit.co/) platform was built, to simplify the way to access and analyze data.

Its aim is to simplify the analysis from mouse/cell cluture samples, which are often a small size. In detailed, the number of sample is limited below 52, due to it hosted on a server with only 4 Cores CPU + 8GB RAM. But if you want to analysis more, please modify code.

Also, due to the resource limitation, password is required for analysis. The password is stored in `Pwd.rds` using a simple string object. To be honest, this approach only deters honest people but not those with ill intentions. However, it does have some effect in conserving lab resources.

# Install

You can simply download all the file, and run using [RStudio](https://posit.co/products/open-source/rstudio/), or host it on [Shiny Server](https://posit.co/products/open-source/shinyserver/).

Before run it, please check `dependency.R` to make sure every package depended on is installed properly.

# Features & Modules Description

## Differential Expression Analysis Module (DEGAnalysis)

-   Basic quality control and filtering (e.g., protein-coding genes and non-protein-coding genes)
-   Differential analysis of bulk RNA-seq (using `DESeq2`, `edgeR`, `limma`)
-   Differential analysis of non-sequencing source data, such as microarrays or other sources (using `limma`)
-   PCA analysis
-   Correlation analysis between samples
-   Sequence trend clustering analysis (`Mfuzz` analysis)
-   Batch effect correction for samples from different sources (using `sva`)

## Sample Scoring Analysis Module (MatAnalysis)

-   Pathway activity inference (`PROGENy` model)
-   Transcription factor activity inference (`CollecTRI` model)
-   `ssGSEA` scoring analysis
-   Immune infiltration scoring (using `xCell`)

## Correlation Analysis Module (CorrAnalysis)

-   Expression correlation analysis between a single gene and all other genes

## Co-expression Module Identification (coExprAnalysis)

-   Identification of gene groups based on non-negative matrix factorization (powered by ultra-fast `RcppML` package)

-   We don't adopt `WGCNA` for its slow calculation, because the limited resource of server.

## Enrichment Analysis Module (TODO)

-   Enrichment analysis via clusterProfiler, supports `ORA`, multi-group `ORA`, and `GSEA`.

-   Results were clustered using `aPEAR` package.

-   Because hosted on server, `KEGG` support is solved by `gson` package.

# Version

Version: 0.1.2 (2024-8-17)
