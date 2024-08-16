library(clusterProfiler)
c2.cp.biocarta<-getGmt("../Msigdb/c2.cp.biocarta.v2024.1.Hs.symbols.gmt")
c2.cp.kegg_legacy<-getGmt("../Msigdb/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")
c2.cp.kegg_medicus<-getGmt("../Msigdb/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt")
c2.cp.pid<-getGmt("../Msigdb/c2.cp.pid.v2024.1.Hs.symbols.gmt")
c2.cp.wikipathways<-getGmt("../Msigdb/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")
h.all<-getGmt("../Msigdb/h.all.v2024.1.Hs.symbols.gmt")
c6.all<-getGmt("../Msigdb/c6.all.v2024.1.Hs.symbols.gmt")
c8.all<-getGmt("~/Downloads/c8.all.v2024.1.Hs.symbols.gmt")


Sys.setenv(http_proxy = "http://127.0.0.1:7890")
Sys.setenv(https_proxy = "http://127.0.0.1:7890")
library(decoupleR)
progenyNet.mouse <- get_progeny(organism = 'mouse', top = 500)
progenyNet.human <- get_progeny(organism = 'human', top = 500)
collectriNet.mouse <- get_collectri(organism = 'mouse', split_complexes=FALSE)
collectriNet.human <- get_collectri(organism = 'human', split_complexes=FALSE)


save(c2.cp.biocarta,c2.cp.kegg_legacy,c2.cp.kegg_medicus,c2.cp.pid,c2.cp.wikipathways,h.all,c6.all,progenyNet.mouse,progenyNet.human,collectriNet.mouse,collectriNet.human,c8.all,file="AnnoPwDBsV2.RData")

