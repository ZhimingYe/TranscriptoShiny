# Transcriptomic analysis on Shiny

## Purpose

I come from a wet lab where most members are not familiar with using R programming for sequencing data analysis. Therefore, this Shiny platform was built to simplify the way everyone can access and analyze data.

## Differential Expression Analysis Module (DEGAnalysis)

-   Basic quality control and filtering (e.g., protein-coding genes and non-protein-coding genes)
-   Differential analysis of bulk RNA-seq (using `DESeq2`, `edgeR`, `limma`)
-   Differential analysis of non-sequencing source data, such as microarrays or other sources (using `limma`)
-   PCA analysis
-   Correlation analysis between samples
-   Sequence trend clustering analysis (`Mfuzz` analysis)
-   Batch effect correction for samples from different sources (using `sva`)

## Sample Scoring Analysis Module (MatAnalysis)

-   Pathway activity inference (PROGENy)
-   Transcription factor activity inference (CollecTRI)
-   `ssGSEA` scoring analysis
-   Immune infiltration scoring (using `xCell`)

## Correlation Analysis Module (CorrAnalysis)

-   Expression correlation analysis between a single gene and all other genes

## Co-expression Module Identification (coExprAnalysis)

-   Identification of gene groups based on non-negative matrix factorization

## Version

Version: 0.1.2 (2024-8-17)
