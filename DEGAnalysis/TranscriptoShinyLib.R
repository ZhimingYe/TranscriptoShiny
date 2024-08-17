read_file. <- function(file_path) {
  file_extension <- tolower(tools::file_ext(file_path))
  if (file_extension == "xlsx") {
    data <- read_excel(file_path)
  } else if (file_extension == "tsv") {
    data <- read.delim(file_path,sep = "\t")
  } else if (file_extension == "csv") {
    data <- read_csv(file_path)
  } else {
    stop("Unsupported file type. Please provide a file with .xlsx, .tsv, or .csv extension.")
  }
  
  return(data)
}


PrefilterDF<-function(ExprTable,Group,SEQtypeDfTypeloggedInput,doBatchremove,NumFilter,FilterPC,Spec){
  SEQtype<-case_when(SEQtypeDfTypeloggedInput%in%c("COUNTRAW","TPMRAW","TPMLOG")~"htseq",
                     SEQtypeDfTypeloggedInput%in%c("LSRAW","LSLOG")~"LS")
  DfType<-case_when(SEQtypeDfTypeloggedInput%in%c("COUNTRAW")~"Count",
                    T~"NOTCount")
  logged<-case_when(SEQtypeDfTypeloggedInput%in%c("COUNTRAW","TPMRAW","LSRAW")~F,
                    T~T)
  
  if(typeof(ExprTable[,1,drop=T])!="character"){
    stop("CSV file format error, can't find gene column.")
  }
  if(sum(grepl("[-]",colnames(ExprTable)))!=0){
    stop("Column names should not contain short dash (-)")
  }
  if((ncol(ExprTable)-1)!=nrow(Group)){
    stop("Expr table not match with the grouping table. [Sample num not matched]")
  }
  if(ncol(ExprTable)>=52){
    stop("Too many samples! Please use computer yourself.")
  }
  if(sum(Group[,1,drop=T]!=colnames(ExprTable)[-1])>0){
    stop("Expr table not match with the grouping table. [Order not matched]")
  }
  showNotification("Start Prefilter")
  if(SEQtype=="htseq"){
    if(nrow(ExprTable)<3000){
      stop("Please input entire dataset")
    }
  }
  if(SEQtype=="LS"){
    if(nrow(ExprTable)<250){
      stop("Please input entire dataset")
    }
  }
  if(sum(sapply(ExprTable, is.numeric)[-1])!=nrow(Group)){
    stop("Expression contains non-numeric")
  }
  TestFrame<-ExprTable[1:100,1]
  ORGan<-"UNKNOWN"
  GEType<-"SYMBOL"
  if(sum(grepl("^ENSMU",TestFrame[,1,drop=T]))>80){
    GEType<-"ENSEMBL"
  }
  if(sum(grepl("^ENSG",TestFrame[,1,drop=T]))>80){
    GEType<-"ENSEMBL"
  }
  if(Spec=="Human"){
    ORGan<-"HUMAN"
  }
  if(Spec=="Mouse"){
    ORGan<-"MOUSE"
  }
  if(Spec=="Others"){
    ORGan<-"HUMAN_OR_UNKNOWN"
  }
  colnames(ExprTable)[1]<-"GeneCol"
  # if(ORGan=="HUMAN_OR_UNKNOWN"){
  #   ExprTable$GeneCol<-str_to_upper(ExprTable$GeneCol)
  # }
  if(GEType=="ENSEMBL"&(sum(grepl("[.]",ExprTable$GeneCol))>10)){
    try({ExprTable$GeneCol<-sapply(strsplit(ExprTable$GeneCol,'[.]'), function(x)x[1])})
  }
  if(sum(duplicated(ExprTable[,1,drop=T]))!=0){
    showNotification("Find duplicated feature, correcting.")
    ExprTable2<-aggregate(ExprTable[,-1],by=list(Gene=ExprTable$GeneCol),FUN=mean,na.action = na.omit)
    ExprTable2<-ExprTable2%>%column_to_rownames(var="Gene")
  }
  else{
    ExprTable2<-ExprTable%>%column_to_rownames(var=colnames(ExprTable)[1])
  }
  if(DfType!="Count"&logged==F){
    ExprTable2<-log2(ExprTable2+1)
  }
  if(SEQtype=="LS"){
    ExprTable2<-normalizeBetweenArrays(ExprTable2)
  }
  if(DfType=="Count"){
    ExprTable2<-round(ExprTable2)
  }
  if(sum(rowSums(ExprTable2>0)>NumFilter)>(nrow(ExprTable2)*0.75)){
    showNotification(paste0("More than half feature low express! will Only keep",(sum(rowSums(ExprTable2>0)>NumFilter))),type = "warning")
    
  }
  ExprTable2<-ExprTable2[rowSums(ExprTable2>0)>NumFilter,]
  
  if(sum(((Group[,2,drop=T])%>%table()%>%as.numeric())>=3)<length((Group[,2,drop=T])%>%table())){
    stop("at least 3 sample each group needed")
  }
  if(doBatchremove){
    showNotification("Removing batch effect...")
    if(sum(((Group[,3,drop=T])%>%table()%>%as.numeric())>=3)<length((Group[,3,drop=T])%>%table())){
      stop("at least 3 sample each batch needed")
    }
    try({
      if(DfType=="Count"){
        ExprTable2<-sva::ComBat_seq(ExprTable2%>%as.matrix(),Group[,3,drop=T]%>%as.factor(),full_mod=F)
      }
      else{
        ExprTable2<-sva::ComBat(ExprTable2%>%as.matrix(),Group[,3,drop=T]%>%as.factor())
      }
    })
  }
  ExprTable2<-as.data.frame(ExprTable2)
  load("~/Documents/4Fun/TranscriptoShiny/DEGAnalysis/GeneDBBiomart.RData")
  ncH<-protein_coding_genesHuman%>%dplyr::filter(gene_biotype!="protein_coding")
  ncM<-protein_coding_genesMouse%>%dplyr::filter(gene_biotype!="protein_coding")
  pcH<-protein_coding_genesHuman%>%dplyr::filter(gene_biotype=="protein_coding")
  pcM<-protein_coding_genesMouse%>%dplyr::filter(gene_biotype=="protein_coding")
  if(FilterPC=="PC"){
    orginrow<-nrow(ExprTable2)
    if(GEType=="SYMBOL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcM$external_gene_name,]
      }
      if(ORGan=="HUMAN"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcH$external_gene_name,]
      }
    }
    if(GEType=="ENSEMBL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcM$ensembl_gene_id,]
      }
      if(ORGan=="HUMAN"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcH$ensembl_gene_id,]
      }
    }
    if(nrow(ExprTable2)<100){
      stop("Unrecognized Organism")
    }
    afterrow<-nrow(ExprTable2)
    showNotification(paste0("Before: ",orginrow,"  After: ",afterrow))
  }
  if(FilterPC=="NC"){
    orginrow<-nrow(ExprTable2)
    if(GEType=="SYMBOL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncM$external_gene_name,]
      }
      if(ORGan=="HUMAN"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncH$external_gene_name,]
      }
    }
    if(GEType=="ENSEMBL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncM$ensembl_gene_id,]
      }
      if(ORGan=="HUMAN"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncH$ensembl_gene_id,]
      }
    }
    if(nrow(ExprTable2)<100){
      stop("Unrecognized Organism")
    }
    afterrow<-nrow(ExprTable2)
    showNotification(paste0("Before: ",orginrow,"  After: ",afterrow))
  }
  # return(ExprTable2)
  # recheck
  if(dim(ExprTable2)[2]==nrow(Group)){
    showNotification("Prefilter Successed.")
    return(ExprTable2)
  }
  else{
    stop("Pre-filter unsuccessful. Please re-check")
  }
}

std_<-function(x){
  data<-x%>%as.matrix()
  for (i in 1:dim(data)[[1]]) {
    data[i,] <- (data[i,] - mean(data[i,], na.rm = TRUE))/sd(data[i,], na.rm = TRUE)
  }
  return(data)
}

Get0ID<-function(x,SpecInput){
  load("./GeneDBBiomart.RData")
  outPut<-x
  TYPEo<-SpecInput
  Flaginternal<-F
  if(sum(grepl("^ENSMU",outPut$GID[1:100]))>80){
    colnames(outPut)[which(colnames(outPut)=="GID")]<-"ensembl_gene_id"
    TYPEo<-"Mouse"
    Flaginternal<-T
  }
  if(sum(grepl("^ENSG",outPut$GID[1:100]))>80){
    colnames(outPut)[which(colnames(outPut)=="GID")]<-"ensembl_gene_id"
    TYPEo<-"Human"
    Flaginternal<-T
  }
  if(!Flaginternal){
    colnames(outPut)[which(colnames(outPut)=="GID")]<-"external_gene_name"
  }
  if(TYPEo=="Human"){
    outPut<-outPut%>%left_join(protein_coding_genesHuman)
  }
  if(TYPEo=="Mouse"){
    outPut<-outPut%>%left_join(protein_coding_genesMouse)
  }
  return(outPut)
  
}

CorrEstimate0 <-
  function(genegsvamatrix,
           TargetSet = "HALLMARK_HYPOXIA",
           UseMethod = "spearman") {
    data <- genegsvamatrix%>%as.data.frame()
    gene_name1 <- c()
    gene_name2 <- c()
    cor_r <- c()
    pvalue <- c()
    i = which(rownames(data) == TargetSet)
    TargetNum<-nrow(data)
    
    for (r in 1:nrow(data)) {
      g1 = rownames(data)[i]
      g2 = rownames(data)[r]
      c_r = cor(as.numeric(data[i, ]), as.numeric(data[r, ]), method = UseMethod)
      p = cor.test(as.numeric(data[i, ]), as.numeric(data[r, ]), method =
                     UseMethod)[[3]]
      gene_name1 = c(gene_name1, g1)
      gene_name2 = c(gene_name2, g2)
      cor_r = c(cor_r, c_r)
      pvalue = c(pvalue, p)
      
      setProgress((r/TargetNum)-0.01)
      
    }
    data_cor <- data.frame(gene_name1, gene_name2, cor_r, pvalue)
    return(data_cor)
  }


CalcMad <- function(mat, num) {
  mads <- apply(mat, 1, mad)
  mat <- mat[rev(order(mads)), ]
  return(mat[1:num, ])
}
PlotPCA<-function(Mat,group,Ellipse=T,counts=F,loged=T,NumOfGenes=5000){
  if(counts){
    Mat<-log2(edgeR::cpm(Mat)+1)
  }
  if((!loged)&(!counts)){
    Mat<-log2(Mat+1)
  }
  if(nrow(Mat)>5000){
    Mat<-CalcMad(Mat,NumOfGenes)
  }
  res.pca2 <- prcomp(Mat%>%t(),scale=T)
  pic1<-ggbiplot.internal(res.pca2,var.axes = F,obs.scale = 0.5,groups = as.factor(group),ellipse = Ellipse,circle = F)+theme_bw()
  showNotification("PCA plot generated")
  return(pic1)
}
PlotCorr<-function(Mat,group,Ellipse=T,counts=F,loged=T,NumOfGenes=5000){
  if(counts){
    Mat<-log2(edgeR::cpm(Mat)+1)
  }
  if((!loged)&(!counts)){
    Mat<-log2(Mat+1)
  }
  if(nrow(Mat)>5000){
    Mat<-CalcMad(Mat,NumOfGenes)
  }
  CORmat<-cor(Mat)
  pic1<-Heatmap(CORmat,name="Corr Coeff.")
  showNotification("Correlation plot generated")
  return(pic1)
}
GetGroups<-function(groups){
  if(length(names(table(groups)))>8)stop("too many groups, please reduce!")
  groups3<-names(table(groups))
  combinations <- expand.grid(groups3, groups3)
  groups2<-names(table(groups))
  all_permutations <- gtools::permutations(n = length(groups2), r = length(groups2), v = groups2)
  all_permutations <- apply(all_permutations, 1, paste, collapse = "-")
  all_permutationsA<-data.frame(Cpr=all_permutations)
  all_permutationsA<-all_permutationsA%>%dplyr::mutate(GroupOrder=1:nrow(all_permutationsA))
  combinations2<-combinations[combinations$Var1!=combinations$Var2,]
  combinations2<-combinations2%>%as.data.frame()
  combinations2<-combinations2%>%dplyr::mutate(Cpr=paste0(Var1,"-",Var2),GroupOrder=1:nrow(combinations2))
  
  return(list(`A`=combinations2,`B`=all_permutationsA))
}
#

GenDEGtable<-function(Expr,GroupInput,CprString){
  CprA<-sapply(strsplit(CprString,"[-]"),function(x)x[1])
  CprB<-sapply(strsplit(CprString,"[-]"),function(x)x[2])
  KeepVec<-GroupInput%in%c(CprA,CprB)
  GroupVec<-GroupInput[KeepVec]
  ExprT<-Expr[,KeepVec]
  
  CprDF<-data.frame(ID=colnames(ExprT),Group=GroupVec)
  showNotification("Trying Perform DESeq2 analysis, Please wait...")
  ddsx<-DESeqDataSetFromMatrix(countData = round(ExprT), colData = CprDF, design = ~Group)
  showNotification("Waiting DESeq2 response... ")
  ddsx <- DESeq(ddsx)
  contrastx <- c("Group", CprA, CprB)
  showNotification("Waiting DESeq2 response... ")
  CPRx <- lfcShrink(ddsx, contrast=contrastx, type="ashr")
  CPRxA <- as.data.frame(CPRx)

  showNotification("Trying Perform edgeR analysis, Please wait...")
  group <- factor(GroupVec)
  y <- DGEList(counts = round(ExprT), group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrast_matrix <- makeContrasts(
    contrasts = CprString,
    levels = design
  )
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  showNotification("Trying Perform limma analysis, Please wait...")
  v <- voom(y, design, plot = F)
  fitV <- lmFit(v, design)
  resultFull1 <- glmQLFTest(fit, contrast=contrast_matrix)
  fit1_voom <- contrasts.fit(fitV, contrast=contrast_matrix)
  fit1v <- eBayes(fit1_voom)
  ResLimma<-topTable(fit1v, coef = paste0(CprA,"-",CprB),number = Inf)
  edgeRresult<-topTags(resultFull1,n = Inf)
  edgeRresult<-as.data.frame(edgeRresult)
  edgeRresult2<-edgeRresult%>%dplyr::select(logFC,PValue,FDR)%>%dplyr::rename(Padj=FDR)%>%dplyr::mutate(Package="edgeR")%>%rownames_to_column(var="GID")
  ResLimma2<-ResLimma%>%dplyr::select(logFC,P.Value,adj.P.Val)%>%dplyr::rename(PValue=P.Value,Padj=adj.P.Val)%>%dplyr::mutate(Package="limma")%>%rownames_to_column(var="GID")
  CPRxA2<-CPRxA%>%dplyr::select(log2FoldChange,pvalue,padj)%>%dplyr::rename(logFC=log2FoldChange,PValue=pvalue,Padj=padj)%>%dplyr::mutate(Package="DESeq2")%>%rownames_to_column(var="GID")
  DEGres<-rbind(CPRxA2,edgeRresult2,ResLimma2)%>%dplyr::arrange(GID)
  return(DEGres)
}



GenDEGtable.Norm<-function(Expr,GroupInput,CprString){
  CprA<-sapply(strsplit(CprString,"[-]"),function(x)x[1])
  CprB<-sapply(strsplit(CprString,"[-]"),function(x)x[2])
  KeepVec<-GroupInput%in%c(CprA,CprB)
  GroupVec<-GroupInput[KeepVec]
  ExprT<-Expr[,KeepVec]
  CprDF<-data.frame(ID=colnames(ExprT),Group=GroupVec)
  showNotification("Trying Perform limma analysis, Please wait...")
  design = model.matrix(~0+GroupVec)
  colnames(design) = gsub("GroupVec","",colnames(design))
  x <- c(CprString)
  contrast =  makeContrasts(contrasts=x,levels=design)
  fit1 <- lmFit(ExprT, design)
  fit2 <- contrasts.fit(fit1,contrasts = contrast)
  fit3 <- eBayes(fit2,trend = T)
  tmpOut <- topTable(fit3,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC")
  return(tmpOut%>%dplyr::select(logFC,P.Value,adj.P.Val)%>%dplyr::rename(PValue=P.Value,Padj=adj.P.Val)%>%dplyr::mutate(Package="limma")%>%rownames_to_column(var="GID"))
}

ConvertEnsembl2Symbol<-function(x,Org){
  df<-x
  if(Org=="Mouse"){
    
    gdf<-bitr(df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
    
    gdf$SYMBOL<-stringr::str_to_upper(gdf$SYMBOL)
    df<-df%>%left_join(gdf)
    return(df)
  }
  else{
    gdf<-bitr(df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
    df<-df%>%left_join(gdf)
    return(df)
  }
}
MouseSymbol2Human<-function(x){
  x$SYMBOL<-stringr::str_to_upper(x$SYMBOL)
  return(x)
}
BuildMultigroupDEGlist<-function(Gene,Group,FromType = "SYMBOL",OrgDB=org.Hs.eg.db,UniversalDB=F){
  Genetable <- data.frame(Gene = Gene, Group = Group)
  Genetable<-Genetable[!is.na(Genetable$Gene),]
  Genetable$testcol<-paste0(Genetable$Gene,"_",Genetable$Group)
  Genetable<-Genetable[!duplicated(Genetable$testcol),]
  ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,
                                         toType = "ENTREZID", OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable)
  if(UniversalDB){
    Genetable$ENTREZID<-Genetable$Gene
  }
  return(Genetable)
}
Ranked.GS<-function(Gene,log2FC,FromType = "SYMBOL",OrgDB=org.Hs.eg.db,UniversalDB=F){
  Genetable<-data.frame(Gene=Gene,log2FC=log2FC)
  Genetable<-Genetable[!is.na(Genetable$Gene),]
  ENTREZIDtable<-clusterProfiler::bitr(Genetable$Gene,fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1]<-"Gene"
  Genetable<-Genetable%>%left_join(ENTREZIDtable)%>%arrange(desc(log2FC))
  GSElist<-as.numeric(Genetable$log2FC)
  if(!UniversalDB){
    names(GSElist)<-Genetable$ENTREZID
  }
  else{
    names(GSElist)<-Genetable$Gene
  }
  GSElist = sort(GSElist, decreasing = TRUE)
  return(GSElist)
}
doGO<-function(MultiGroup="SG",GS,ont="BP",PVal=0.01,QVal=0.05,FromType = "SYMBOL",OrgDB=org.Hs.eg.db){
  if(MultiGroup=="SG"){
    GeneList<-GS%>%unique()%>%bitr(fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
    EnrichGO<-enrichGO(GeneList$ENTREZID,ont = ont,OrgDb = OrgDB,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="MG"){
    EnrichGO<-compareCluster(ENTREZID~Group, data=GS, fun="enrichGO",ont = ont,OrgDb = OrgDB,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="GSEA"){
    EnrichGO<-gseGO(GS,pvalueCutoff = PVal,ont = ont,OrgDb = OrgDB)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  return(EnrichGO)
}
doKEGG<-function(MultiGroup="SG",GS,PVal=0.01,QVal=0.05,FromType = "SYMBOL",Organism="hsa",OrgDB=org.Hs.eg.db,GSONFILE){
  library(gson)
  KEGGgson<-GSONFILE
  if(MultiGroup=="SG"){
    GeneList<-GS%>%unique()%>%bitr(fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
    EnrichKEGG<-enricher(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal,gson = KEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="MG"){
    EnrichKEGG<-compareCluster(ENTREZID~Group, data=GS, fun="enricher",pvalueCutoff = PVal,qvalueCutoff = QVal,gson = KEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="GSEA"){
    EnrichKEGG<-GSEA(GS,pvalueCutoff = PVal,gson=KEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  return(EnrichKEGG)
}
doMKEGG<-function(MultiGroup="SG",GS,PVal=0.01,QVal=0.05,FromType = "SYMBOL",Organism="hsa",OrgDB=org.Hs.eg.db,GSONFILE){
  library(gson)
  MKEGGgson<-GSONFILE
  if(MultiGroup=="SG"){
    GeneList<-GS%>%unique()%>%bitr(fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
    EnrichKEGG<-enricher(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal,gson = MKEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="MG"){
    EnrichKEGG<-compareCluster(ENTREZID~Group, data=GS, fun="enricher",pvalueCutoff = PVal,qvalueCutoff = QVal,gson = MKEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="GSEA"){
    EnrichKEGG<-GSEA(GS,pvalueCutoff = PVal,gson=MKEGGgson)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  return(EnrichKEGG)
}
setReadable2<-function(x,OrgDb,keyType){
  x@readable<-F
  return(clusterProfiler::setReadable(x,OrgDb,keyType))
}
doRA<-function(MultiGroup,GS,PVal=0.01,QVal=0.05,FromType = "SYMBOL",Organism="human",OrgDB=org.Hs.eg.db){
  library(ReactomePA)
  if(MultiGroup=="SG"){
    GeneList<-GS%>%unique()%>%bitr(fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
    EnrichRA<-enrichPathway(GeneList$ENTREZID,pvalueCutoff = PVal,qvalueCutoff = QVal,organism=Organism)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="MG"){
    EnrichRA<-compareCluster(ENTREZID~Group, data=GS, fun="enrichPathway",pvalueCutoff = PVal,qvalueCutoff = QVal,organism=Organism)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="GSEA"){
    EnrichRA<-gsePathway(GS,pvalueCutoff = PVal,organism=Organism)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  return(EnrichRA)
}
doDO<-function(MultiGroup,GS,ont="DO",PVal=0.01,QVal=0.05,FromType = "SYMBOL",OrgDB=org.Hs.eg.db){
  if(MultiGroup=="SG"){
    GeneList<-GS%>%unique()%>%bitr(fromType = FromType,toType = "ENTREZID",OrgDb = OrgDB)
    EnrichGO<-enrichDO(GeneList$ENTREZID,ont = ont,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="MG"){
    EnrichGO<-compareCluster(ENTREZID~Group, data=GS, fun="enrichDO",ont = ont,pvalueCutoff = PVal,qvalueCutoff = QVal)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  if(MultiGroup=="GSEA"){
    EnrichGO<-gseDO(GS,pvalueCutoff = PVal)%>%setReadable2(OrgDb=OrgDB,keyType = "ENTREZID")
  }
  return(EnrichGO)
}
doUniversal<-function(MultiGroup="SG",GS,PVal=0.01,QVal=0.05,FromType = "SYMBOL",DB,OrgDB=org.Hs.eg.db){
  if(MultiGroup=="SG"){
    EnrichKEGG<-enricher(GS,pvalueCutoff = PVal,qvalueCutoff = QVal,TERM2GENE = DB)
  }
  if(MultiGroup=="MG"){
    EnrichKEGG<-compareCluster(ENTREZID~Group, data=GS, fun="enricher",pvalueCutoff = PVal,qvalueCutoff = QVal,TERM2GENE = DB)
  }
  if(MultiGroup=="GSEA"){
    EnrichKEGG<-GSEA(GS,pvalueCutoff = PVal,TERM2GENE = DB)
  }
  return(EnrichKEGG)
}
# doRA(MultiGroup = "MG",GS = GSt)->kk


ggbiplot.internal <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                              obs.scale = 1 - scale, var.scale = scale,
                              groups = NULL, ellipse = FALSE, ellipse.prob = 0.68,
                              labels = NULL, labels.size = 3, alpha = 1,
                              var.axes = TRUE,
                              circle = FALSE, circle.prob = 0.69,
                              varname.size = 3, varname.adjust = 1.5,
                              varname.abbrev = FALSE, ...)
{

  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs,
                       sprintf('(%0.1f%% explained var.)',
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) +
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle)
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'),
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')),
                   color = muted('red'))
  }
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  if(var.axes) {
    g <- g +
      geom_text(data = df.v,
                aes(label = varname, x = xvar, y = yvar,
                    angle = angle, hjust = hjust),
                color = 'darkred', size = varname.size)
  }
  return(g)
}



# 

# library(biomaRt)
# Sys.setenv(http_proxy="http://127.0.0.1:7890")
# Sys.setenv(https_proxy="http://127.0.0.1:7890")
# 
# ensemblMouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# gene_LENGTHinfoMouse <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", 'start_position', 'end_position'),mart = ensemblMouse)
# protein_coding_genesMouse <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name","hgnc_symbol", "gene_biotype","description"),
#   mart = ensemblMouse
# )
# 
# ensemblHuman <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_LENGTHinfoHuman <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", 'start_position', 'end_position'),mart = ensemblHuman)
# protein_coding_genesHuman <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name","hgnc_symbol", "gene_biotype","description"),
#   mart = ensemblHuman
# )
# save(gene_LENGTHinfoMouse,protein_coding_genesMouse,gene_LENGTHinfoHuman,protein_coding_genesHuman,file="GeneDBBiomart.RData")
