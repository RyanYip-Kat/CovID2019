library(Seurat)
library(dplyr)
library(stringr)


t<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/VHK_20200307_9/monocle3/model/20200309_seurat.rds")
c<-readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/all_gene/PCV_CAT_8/monocle3/model/20200309_seurat.rds")

Idents(t)=t$new_umap_clusters
Idents(c)=c$new_umap_clusters
t<-subset(t,idents=c(1,6))
c<-subset(c,idents=c(1,7))
seurat<-RunCCA(t,c,renormalize=TRUE,rescale=TRUE)

idents<-ifelse(str_detect(seurat$time,"UVR|UVQ"),"VKH","PCVCAT")
Idents(seurat)<-idents
markers<-FindMarkers(seurat,ident.1="PCVCAT",ident.2="VKH")
cluster<-with(markers,ifelse(avg_logFC>0,"PCVCAT","VKH"))
markers$cluster=cluster
write.table(markers,"vhk_pcvcat_markers.csv",sep=",",quote=F)



