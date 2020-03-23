library(Seurat)
library(monocle3)
library(dplyr)
library(argparse)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="")
parser$add_argument("--cds",
                      type="character",
                      default="",
                      help="")

parser$add_argument("-k",
                      type="integer",
                      default="20",
                      help="")

parser$add_argument("--outdir",
		     type="character",
		     default="")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
cds<-readRDS(args$cds)

Idents(seurat)<-seurat$UMAP_clusters
#seurat<-subset(seurat,idents=c(1, 2, 4, 5, 8, 9,15,12,19,36)) # VKH
seurat<-subset(seurat,idents=c(2, 3, 5, 7, 8, 17, 18, 22, 24))
#seurat<-FindVariableFeatures(pbmc,nfeatures=2000)
cells<-Cells(seurat)
cds<-cds[,cells]
cds<-cluster_cells(cds,
                   reduction_method="UMAP",
                   k=as.integer(args$k),
                   cluster_method="leiden",
                   partition_qval=0.01)


UMAP_clusters<-clusters(cds,reduction_method="UMAP")
UMAP_clusters<-as.data.frame(UMAP_clusters)
colnames(UMAP_clusters)<-"new_umap_clusters"

seurat<-AddMetaData(seurat,metadata=UMAP_clusters)
saveRDS(seurat,file.path(args$outdir,"20200309_seurat.rds"))
jpeg(file.path(args$outdir,"new_umap_clusters.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="umap",group.by="new_umap_clusters",label=TRUE,label.size=7.5)
dev.off()

Idents(seurat)<-seurat$new_umap_clusters
markers <- FindAllMarkers(seurat, only.pos = FALSE,
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )
topn=markers%>%group_by(cluster)%>%top_n(70,wt=avg_logFC)
saveRDS(markers,file.path(args$outdir,"new_markers.rds"))
write.table(topn,file.path(args$outdir,"new_top70_markers.csv"),quote=F,row.names=F,sep=",")
