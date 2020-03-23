library(Seurat)
library(dplyr)
library(stringr)


t<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/VHK_20200307_9/monocle3/model/seurat.rds")
c<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/PCV_CAT_8/monocle3/model/seurat.rds")

Idents(t)=t$UMAP_clusters
Idents(c)=c$UMAP_clusters
t<-subset(t,idents=c(1, 2, 4, 5, 8, 9,15,12,36))
c<-subset(c,idents=c(2, 3, 5, 7, 8, 17, 18, 22, 24))
seurat<-RunCCA(t,c,renormalize=TRUE,rescale=TRUE)

idents<-ifelse(str_detect(seurat$time,"UVR|UVQ"),"VKH","PCVCAT")
Idents(seurat)<-idents
markers<-FindMarkers(seurat,ident.1="PCVCAT",ident.2="VKH")
cluster<-with(markers,ifelse(avg_logFC>0,"PCVCAT","VKH"))
markers$cluster=cluster
write.table(markers,"Su/20200309/vhk_pcvcat_markers.csv",sep=",",quote=F)



