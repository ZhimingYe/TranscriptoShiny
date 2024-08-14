

tst<-function(){
  if(typeof(qzhFPKM[,1,drop=T])!="character"){
    stop("CSV file format error, can't find gene column.")
  }
  if((ncol(qzhFPKM)-1)!=nrow(Group)){
    stop("Expr table not match with the grouping table. [Sample num not matched]")
  }
  if(Group[,1,drop=T]!=colnames(qzhFPKM)[-1]){
    stop("Expr table not match with the grouping table. [Order not matched]")
  }
  if(SEQtype=="htseq"){
    if(nrow(qzhFPKM)<3000){
      stop("Please input entire dataset")
    }
  }
  if(SEQtype=="LS"){
    if(nrow(qzhFPKM)<250){
      stop("Please input entire dataset")
    }
  }
  if(sum(sapply(qzhFPKM, is.numeric)[-1])!=nrow(Group)){
    stop("Expression contains non-numeric")
  }
  TestFrame<-qzhFPKM[1:10,1]
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
    ORGan<-"HUMAN_OR_UNKNOWN"
    GEType<-"SYMBOL"
  }
  colnames(qzhFPKM)[1]<-GeneCol
  if(ORGan=="HUMAN_OR_UNKNOWN"){
    qzhFPKM$GeneCol<-str_to_upper(qzhFPKM$GeneCol)
  }
  if(sum(duplicated(qzhFPKM[,1,drop=T]))!=0){
    qzhFPKM2<-aggregate(qzhFPKM[,-1],by=list(Gene=qzhFPKM$GeneCol),FUN=mean,na.action = na.omit)
    qzhFPKM2<-qzhFPKM2%>%column_to_rownames(var="Gene")
  }
  else{
    qzhFPKM2<-qzhFPKM%>%column_to_rownames(var=colnames(qzhFPKM)[1])
  }
  if(DfType!="Count"&logged==F){
    qzhFPKM2<-log2(qzhFPKM2+1)
  }
  if(SEQtype=="LS"){
    qzhFPKM2<-normalizeBetweenArrays(qzhFPKM2)
  }
  if(DfType=="Count"){
    qzhFPKM2<-round(qzhFPKM2)
  }
  qzhFPKM2<-[rowSums(qzhFPKM2>0)>NumFilter,]
  if(sum(((Group[,2,drop=T])%>%table()%>%as.numeric())>=3)>=2){
    stop("at least 3 sample each group needed")
  }
  if(doBatchremove){
    if(sum(((Group[,3,drop=T])%>%table()%>%as.numeric())>=3)>=2){
      stop("at least 3 sample each batch needed")
    }
    try({
      if(DfType=="Count"){
        qzhFPKM2<-sva::ComBat_seq(qzhFPKM2,Group[,3,drop=T])
      }
      else{
        qzhFPKM2<-sva::ComBat(qzhFPKM2,Group[,3,drop=T])
      }
    })
  }
  
  # recheck
  # FILTERING
    # TABLE OBS
    # GET PC
    #GET COMBAT
  
}

tst()

Test<-qzhFPKM[1:10,1]

sum((c("A","A","A","A","B","B","B")%>%table()%>%as.numeric())>=3)>=2

%>%table()%>%as.numeric()


sum(duplicated(qzhFPKM[,1,drop=T]))


apply(qzhFPKM,2,is.numeric)


strings <- c("Hello123", "World45", "TestCase1", "example12", "Invali1dstring", "R2D2")

is_proper_case_with_digits <- grepl("^[A-Z][a-z0-9]*$", strings)

is_proper_case_with_digits

# library(biomaRt)
# Sys.setenv(http_proxy="http://127.0.0.1:7890")
# Sys.setenv(https_proxy="http://127.0.0.1:7890")
# 
# ensemblMouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# gene_LENGTHinfoMouse <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", 'start_position', 'end_position'),mart = ensemblMouse)
# protein_coding_genesMouse <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name","hgnc_symbol", "gene_biotype"),
#   mart = ensemblMouse
# )
# 
# ensemblHuman <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_LENGTHinfoHuman <- getBM(attributes = c('ensembl_gene_id',"hgnc_symbol", 'start_position', 'end_position'),mart = ensemblHuman)
# protein_coding_genesHuman <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name","hgnc_symbol", "gene_biotype"),
#   mart = ensemblHuman
# )
# save(gene_LENGTHinfoMouse,protein_coding_genesMouse,gene_LENGTHinfoHuman,protein_coding_genesHuman,file="GeneDBBiomart.RData")
