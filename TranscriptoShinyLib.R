

tst<-function(ExprTable,Group,SEQtype,DfType,logged,doBatchremove,FilterPC){
  if(typeof(ExprTable[,1,drop=T])!="character"){
    stop("CSV file format error, can't find gene column.")
  }
  if((ncol(ExprTable)-1)!=nrow(Group)){
    stop("Expr table not match with the grouping table. [Sample num not matched]")
  }
  if(Group[,1,drop=T]!=colnames(ExprTable)[-1]){
    stop("Expr table not match with the grouping table. [Order not matched]")
  }
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
  TestFrame<-ExprTable[1:10,1]
  ORGan<-"UNKNOWN"
  GEType<-"SYMBOL"
  if(sum(grepl("^ENSMU",TestFrame[,1,drop=T]))!=10){
    GEType<-"ENSEMBL"
    ORGan<-"MOUSE"
  }
  if(sum(grepl("^ENSG",TestFrame[,1,drop=T]))!=10){
    GEType<-"ENSEMBL"
    ORGan<-"HUMAN"
  }
  if((sum(grepl("^[A-Z][a-z0-9]*$",TestFrame[,1,drop=T]))==10)&GEType=="SYMBOL"){
    ORGan<-"MOUSE"
    GEType<-"SYMBOL"
  }
  else{
    ORGan<-"HUMAN"
    GEType<-"SYMBOL"
  }
  colnames(ExprTable)[1]<-GeneCol
  if(ORGan=="HUMAN_OR_UNKNOWN"){
    ExprTable$GeneCol<-str_to_upper(ExprTable$GeneCol)
  }
  if(GEType=="ENSEMBL"&(sum(grepl("[.]",ExprTable$GeneCol))>10)){
    try({ExprTable$GeneCol<-sapply(strsplit(ExprTable$GeneCol,'[.]'), function(x)x[1])})
  }
  if(sum(duplicated(ExprTable[,1,drop=T]))!=0){
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
  ExprTable2<-[rowSums(ExprTable2>0)>NumFilter,]
  if(sum(((Group[,2,drop=T])%>%table()%>%as.numeric())>=3)>=2){
    stop("at least 3 sample each group needed")
  }
  if(doBatchremove){
    if(sum(((Group[,3,drop=T])%>%table()%>%as.numeric())>=3)>=2){
      stop("at least 3 sample each batch needed")
    }
    try({
      if(DfType=="Count"){
        ExprTable2<-sva::ComBat_seq(ExprTable2,Group[,3,drop=T])
      }
      else{
        ExprTable2<-sva::ComBat(ExprTable2,Group[,3,drop=T])
      }
    })
  }
  load("~/Documents/4Fun/TranscriptoShiny/GeneDBBiomart.RData")
  ncH<-protein_coding_genesHuman%>%dplyr::filter(gene_biotype!="protein_coding")
  ncM<-protein_coding_genesMouse%>%dplyr::filter(gene_biotype!="protein_coding")
  pcH<-protein_coding_genesHuman%>%dplyr::filter(gene_biotype=="protein_coding")
  pcM<-protein_coding_genesMouse%>%dplyr::filter(gene_biotype=="protein_coding")
  if(FilterPC=="PC"){
    if(GEType=="SYMBOL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcM$external_gene_name,]
      }
      else{
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcH$external_gene_name,]
      }
    }
    if(GEType=="ENSEMBL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcM$ensembl_gene_id,]
      }
      else{
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%pcH$ensembl_gene_id,]
      }
    }
    if(nrow(ExprTable2)<100){
      stop("Unrecognized Organism")
    }
  }
  if(FilterPC=="NC"){
    if(GEType=="SYMBOL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncM$external_gene_name,]
      }
      else{
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncH$external_gene_name,]
      }
    }
    if(GEType=="ENSEMBL"){
      if(ORGan=="MOUSE"){
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncM$ensembl_gene_id,]
      }
      else{
        ExprTable2<-ExprTable2[rownames(ExprTable2)%in%ncH$ensembl_gene_id,]
      }
    }
    if(nrow(ExprTable2)<100){
      stop("Unrecognized Organism")
    }
  }
  # recheck
  if(dim(ExprTable2)[2]==nrow(Group)){
    return(ExprTable2)
  }
  else{
    stop("Pre-filter unsuccessful. Please re-check")
  }
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
